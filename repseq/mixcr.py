from .common_functions import print_progress_bar, filter_by_functionality
import pandas as pd
import numpy as np
import os
import shlex
import subprocess
from time import sleep
from .slurm import run_slurm_command_from_jupyter
from .io import open_json_report, read_json_report, read_clonoset
from .clonosets import find_all_exported_clonosets_in_folder
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import seaborn as sns
import re
import warnings

JOB_TABLE_COLUMNS = [
    "jobname",
    "sample_id",
    "backend",
    "job_id",
    "status",
    "returncode",
    "cwd",
    "log_filename",
    "command",
]

DEFAULT_SAMPLE_BARCODE_TAG_PATTERN = (
    r"^N{0:2}tggtatcaacgcagagt(SMPL:N{5})(UMI:N{14})N{1}gctN{16}(R1:*)\^N{20}(R2:*)"
)


def _validate_backend(backend):
    if backend not in ["local", "slurm"]:
        raise ValueError("backend must be either 'local' or 'slurm'")


def _normalize_memory(memory, min_memory=16, max_memory=1500):
    if not isinstance(memory, int):
        raise TypeError("memory parameter must be an integer")
    if memory < min_memory:
        print(f"{memory} < than limit ({min_memory}), using {min_memory} GB")
        return min_memory
    if memory > max_memory:
        print(f"{memory} > than limit ({max_memory}), using {max_memory} GB")
        return max_memory
    return memory


def _strip_mixcr_template_placeholders(command_template):
    remove_list = {"mixcr", "r1", "r2", "output_prefix"}
    return [token for token in shlex.split(command_template) if token not in remove_list]


def _mixcr_analyze_command(mixcr_path, memory, command_template_parts, r1, r2,
                           output_prefix, tag_pattern=None, extra_options=None):
    command_parts = [mixcr_path, f"-Xmx{memory}g", *command_template_parts]
    if tag_pattern is not None:
        command_parts.extend(["--tag-pattern", str(tag_pattern)])
    if extra_options is not None:
        command_parts.extend(extra_options)
    command_parts.extend([str(r1), str(r2), str(output_prefix)])
    return shlex.join(command_parts)


def _sample_barcoded_groups(sample_df):
    groups = {}
    grouped_positions = set()

    for _, group in sample_df.groupby(["R1", "R2"], sort=False, dropna=False):
        if group["sample_id"].nunique(dropna=False) <= 1:
            continue

        if "mix_id" not in sample_df.columns:
            raise ValueError("sample_df must contain a 'mix_id' column for sample-barcoded files")
        if group["mix_id"].isna().any() or group["mix_id"].nunique(dropna=False) != 1:
            raise ValueError("All rows with the same R1 and R2 files must have the same non-empty mix_id")
        mix_id = group["mix_id"].iloc[0]
        if isinstance(mix_id, str) and not mix_id.strip():
            raise ValueError("mix_id must be non-empty for sample-barcoded files")

        if "SMPL" not in sample_df.columns:
            raise ValueError("sample_df must contain an 'SMPL' column for sample-barcoded files")
        sample_barcodes = group["SMPL"].tolist()
        if not all(isinstance(barcode, str) for barcode in sample_barcodes):
            raise ValueError(f"All SMPL values for mix_id '{mix_id}' must be strings")
        if len(set(sample_barcodes)) != len(sample_barcodes):
            raise ValueError(f"SMPL values for mix_id '{mix_id}' must be distinct")
        if len({len(barcode) for barcode in sample_barcodes}) != 1:
            raise ValueError(f"All SMPL values for mix_id '{mix_id}' must have the same length")
        if not all(re.fullmatch(r"[ATGC]+", barcode) for barcode in sample_barcodes):
            raise ValueError(f"SMPL values for mix_id '{mix_id}' may contain only ATGC letters")

        tag_patterns = []
        if "miNNNPattern" in sample_df.columns:
            for tag_pattern in group["miNNNPattern"]:
                if pd.isna(tag_pattern) or (isinstance(tag_pattern, str) and not tag_pattern.strip()):
                    continue
                if not isinstance(tag_pattern, str):
                    raise ValueError(f"miNNNPattern values for mix_id '{mix_id}' must be strings or empty")
                tag_patterns.append(tag_pattern)
        unique_tag_patterns = list(dict.fromkeys(tag_patterns))
        if len(unique_tag_patterns) > 1:
            raise ValueError(f"All non-empty miNNNPattern values for mix_id '{mix_id}' must be the same")
        tag_pattern = unique_tag_patterns[0] if unique_tag_patterns else DEFAULT_SAMPLE_BARCODE_TAG_PATTERN

        positions = group.index.tolist()
        groups[positions[0]] = {
            "rows": group,
            "mix_id": mix_id,
            "tag_pattern": tag_pattern,
        }
        grouped_positions.update(positions)

    return groups, grouped_positions


def _job_table(jobs, backend):
    rows = []
    for job in jobs:
        rows.append({
            "jobname": job["jobname"],
            "sample_id": job.get("sample_id", ""),
            "backend": backend,
            "job_id": "",
            "status": "pending",
            "returncode": "",
            "cwd": job["cwd"],
            "log_filename": job["log_filename"],
            "command": job["command"],
        })
    return pd.DataFrame(rows, columns=JOB_TABLE_COLUMNS).astype(object)


def _write_batch_table(filename, table, program_name=None):
    if os.path.exists(filename):
        os.remove(filename)
    with open(filename, "w") as f:
        if program_name is not None:
            f.write(f"# {program_name}\n")
        table.to_csv(f, sep="\t", index=False)


def _result_table_filename(batch_filename):
    logs_folder = os.path.join(os.path.dirname(batch_filename), "logs")
    return os.path.join(logs_folder, f"{os.path.splitext(os.path.basename(batch_filename))[0]}_jobs.csv")


def _uses_result_table(batch_filename):
    return os.path.basename(batch_filename) != "mixcr_reports_slurm_batch.log"


def _write_result_table(batch_filename, table):
    if not _uses_result_table(batch_filename):
        return None
    table_filename = _result_table_filename(batch_filename)
    os.makedirs(os.path.dirname(table_filename), exist_ok=True)
    existing_table = _read_result_table(batch_filename)
    merged_table = _merge_result_table(existing_table, table)
    merged_table.to_csv(table_filename, index=False)
    return table_filename


def _read_result_table(batch_filename):
    table_filename = _result_table_filename(batch_filename)
    if not os.path.exists(table_filename):
        return pd.DataFrame(columns=JOB_TABLE_COLUMNS)
    table = pd.read_csv(table_filename, keep_default_na=False)
    for column in JOB_TABLE_COLUMNS:
        if column not in table.columns:
            table[column] = ""
    return table[JOB_TABLE_COLUMNS]


def _merge_result_table(existing_table, update_table):
    if len(existing_table) == 0:
        return update_table[JOB_TABLE_COLUMNS].copy()
    rerun_jobnames = set(update_table["jobname"])
    kept_table = existing_table.loc[~existing_table["jobname"].isin(rerun_jobnames)]
    return pd.concat([kept_table, update_table[JOB_TABLE_COLUMNS]], ignore_index=True)


def _read_batch_table(filename):
    with open(filename, "r") as f:
        first_line = f.readline().rstrip("\n")
    program_name = first_line[2:] if first_line.startswith("# ") else "MiXCR Batch"
    table = pd.read_csv(filename, sep="\t", comment="#", keep_default_na=False)
    return program_name, table


def _resolve_batch_filename(path, default_filename="mixcr_analyze_batch.log"):
    if os.path.isdir(path):
        filename = os.path.join(path, default_filename)
        if not os.path.exists(filename):
            batch_files = [
                os.path.join(path, f)
                for f in os.listdir(path)
                if f.startswith("mixcr_") and f.endswith("_batch.log")
            ]
            if len(batch_files) == 1:
                filename = batch_files[0]
        return filename
    return path


def _status_counts(table):
    finished = int((table["status"] == "finished").sum())
    failed = int(table["status"].isin(["failed", "submit_failed"]).sum())
    total = len(table)
    return finished, failed, total


def _terminal_statuses():
    return ["finished", "failed", "submit_failed"]


def _map_slurm_state(state):
    state = str(state).upper().split()[0]
    if state in ["PENDING", "CONFIGURING", "REQUEUED", "RESIZING", "SUSPENDED"]:
        return "pending"
    if state in ["RUNNING", "COMPLETING", "STAGE_OUT", "SIGNALING", "STOPPED"]:
        return "running"
    if state == "COMPLETED":
        return "finished"
    if state in [
        "BOOT_FAIL",
        "CANCELLED",
        "DEADLINE",
        "FAILED",
        "NODE_FAIL",
        "OUT_OF_MEMORY",
        "PREEMPTED",
        "REVOKED",
        "SPECIAL_EXIT",
        "TIMEOUT",
    ]:
        return "failed"
    return str(state).lower()


def _query_squeue(job_ids):
    if len(job_ids) == 0:
        return {}
    command = ["squeue", "-h", "-j", ",".join(job_ids), "-o", "%A\t%T"]
    try:
        process = subprocess.run(command, capture_output=True, text=True, check=False)
    except FileNotFoundError:
        return {}
    if process.returncode != 0:
        return {}

    statuses = {}
    for line in process.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) < 2:
            continue
        job_id, state = parts[0].strip(), parts[1].strip()
        statuses[job_id] = {"status": _map_slurm_state(state), "returncode": ""}
    return statuses


def _query_sacct(job_ids):
    if len(job_ids) == 0:
        return {}

    statuses = {}
    for requested_job_id in job_ids:
        for extra_args in [["-X"], []]:
            command = [
                "sacct",
                "-n",
                "-P",
                *extra_args,
                "-j",
                str(requested_job_id),
                "-o",
                "JobID,State,ExitCode",
            ]
            try:
                process = subprocess.run(command, capture_output=True, text=True, check=False)
            except FileNotFoundError:
                return statuses
            if process.returncode != 0:
                continue
            for line in process.stdout.splitlines():
                parts = line.split("|")
                if len(parts) < 3:
                    continue
                raw_job_id, state, exit_code = [p.strip() for p in parts[:3]]
                job_id = raw_job_id.split(".")[0]
                if job_id != str(requested_job_id):
                    continue
                status = _map_slurm_state(state)
                if status not in ["finished", "failed"]:
                    continue
                returncode = exit_code.split(":")[0] if exit_code else ""
                statuses[job_id] = {"status": status, "returncode": returncode}
            if str(requested_job_id) in statuses:
                break
    return statuses


def _query_scontrol(job_ids):
    statuses = {}
    for job_id in job_ids:
        command = ["scontrol", "show", "job", str(job_id)]
        try:
            process = subprocess.run(command, capture_output=True, text=True, check=False)
        except FileNotFoundError:
            return statuses
        if process.returncode != 0:
            continue
        state_match = re.search(r"\bJobState=([A-Z_]+)", process.stdout)
        exit_match = re.search(r"\bExitCode=([0-9]+):[0-9]+", process.stdout)
        if state_match is None:
            continue
        status = _map_slurm_state(state_match.group(1))
        returncode = exit_match.group(1) if exit_match else ""
        statuses[str(job_id)] = {"status": status, "returncode": returncode}
    return statuses


def _sync_slurm_status(filename, program_name, table):
    changed = False
    active_mask = (
        (table["backend"] == "slurm")
        & (table["job_id"].astype(str) != "")
        & (~table["status"].isin(_terminal_statuses()))
    )
    job_ids = list(table.loc[active_mask, "job_id"].astype(str))
    if len(job_ids) == 0:
        return table

    statuses = _query_squeue(job_ids)
    missing_job_ids = [job_id for job_id in job_ids if job_id not in statuses]
    statuses.update(_query_sacct(missing_job_ids))
    missing_job_ids = [job_id for job_id in job_ids if job_id not in statuses]
    statuses.update(_query_scontrol(missing_job_ids))
    missing_job_ids = [job_id for job_id in job_ids if job_id not in statuses]
    statuses.update(_query_slurm_logs(table, missing_job_ids))

    for idx, row in table.loc[active_mask].iterrows():
        job_status = statuses.get(str(row["job_id"]))
        if job_status is None:
            continue
        for key, value in job_status.items():
            if str(table.loc[idx, key]) != str(value):
                table.loc[idx, key] = value
                changed = True
    if changed:
        _write_batch_table(filename, table, program_name=program_name)
        if _uses_result_table(filename):
            _write_result_table(filename, table)
    return table


def _query_slurm_logs(table, job_ids):
    statuses = {}
    for _, row in table.iterrows():
        job_id = str(row.get("job_id", ""))
        if job_id not in job_ids:
            continue
        log_filename = row.get("log_filename", "")
        if not isinstance(log_filename, str) or not os.path.exists(log_filename):
            continue
        try:
            with open(log_filename, "r", errors="replace") as f:
                text = f.read()
        except OSError:
            continue
        if "__REPSEQ_STATUS__:finished" in text:
            statuses[job_id] = {"status": "finished", "returncode": "0"}
        else:
            failed_match = re.search(r"__REPSEQ_STATUS__:failed:([0-9]+)", text)
            if failed_match is not None:
                statuses[job_id] = {"status": "failed", "returncode": failed_match.group(1)}
            elif "__REPSEQ_STATUS__:failed" in text:
                statuses[job_id] = {"status": "failed", "returncode": "1"}
    return statuses


def _print_failed_jobs(table):
    failed_table = table.loc[table["status"].isin(["failed", "submit_failed"])]
    if len(failed_table) == 0:
        return
    print("\nFailed jobs:")
    for _, row in failed_table.iterrows():
        print(f"{row['jobname']}\t{row['job_id']}")


def check_batch_progress(path, loop=False, default_filename="mixcr_analyze_batch.log"):
    filename = _resolve_batch_filename(path, default_filename=default_filename)

    while True:
        program_name, table = _read_batch_table(filename)
        table = _sync_slurm_status(filename, program_name, table)
        finished, failed, tasks = _status_counts(table)
        if not loop:
            print(table[["jobname", "sample_id", "backend", "job_id", "status", "returncode", "log_filename"]].to_string(index=False))
        print_progress_bar(finished + failed, tasks, program_name=program_name, object_name="task(s)")
        if not loop or finished + failed == tasks:
            break
        sleep(1)
    _print_failed_jobs(table)


def _update_job_status(batch_filename, program_name, table, jobname, **updates):
    idx = table.index[table["jobname"] == jobname]
    if len(idx) != 1:
        raise ValueError(f"Could not find job '{jobname}' in batch table")
    for key, value in updates.items():
        table.loc[idx[0], key] = value
    _write_batch_table(batch_filename, table, program_name=program_name)


def _run_local_command(job, program_name, batch_filename, table):
    jobname = job["jobname"]
    command = job["command"]
    cwd = job["cwd"]
    log_filename = job["log_filename"]
    _update_job_status(batch_filename, program_name, table, jobname, status="running")
    with open(log_filename, "w") as log:
        process = subprocess.run(
            command,
            shell=True,
            cwd=cwd,
            stdout=log,
            stderr=subprocess.STDOUT,
            text=True,
        )
    status = "finished" if process.returncode == 0 else "failed"
    _update_job_status(
        batch_filename,
        program_name,
        table,
        jobname,
        status=status,
        returncode=process.returncode,
    )


def _parse_slurm_job_id(stdout):
    text = stdout.decode(errors="replace") if isinstance(stdout, bytes) else str(stdout)
    match = re.search(r"Submitted batch job\s+([0-9]+)", text)
    return match.group(1) if match else ""


def _submit_slurm_command(job, cpus, time_estimate, memory, constraint=None):
    command = job["command"]
    jobname = job["jobname"]
    cwd = job["cwd"]
    log_filename = job["log_filename"]
    slurm_command = command
    if cwd is not None:
        slurm_command = f"cd {shlex.quote(cwd)} && {command}"
    slurm_command = (
        f"({slurm_command}); "
        f"status=$?; "
        f'if [ "$status" -eq 0 ]; then echo "__REPSEQ_STATUS__:finished"; '
        f'else echo "__REPSEQ_STATUS__:failed:$status"; fi; '
        f'exit "$status"'
    )
    stdout, stderr = run_slurm_command_from_jupyter(
        slurm_command,
        jobname,
        cpus,
        time_estimate,
        memory,
        log_filename=log_filename,
        verbose=False,
        constraint=constraint,
    )
    return _parse_slurm_job_id(stdout), stdout, stderr


def _save_result_table(table, table_filename):
    if table_filename is None:
        return
    print(f"Logs folder: {os.path.dirname(table_filename)}")
    print(f"Job table: {table_filename}")


def _run_mixcr_jobs(jobs, program_name, batch_filename, backend="local",
                    cpus=40, time_estimate=1.5, memory=32, save_result_table=True,
                    constraint=None):
    _validate_backend(backend)
    table = _job_table(jobs, backend)
    _write_batch_table(batch_filename, table, program_name=program_name)
    if save_result_table:
        _write_result_table(batch_filename, table)
    else:
        table_filename = _result_table_filename(batch_filename)
        if os.path.exists(table_filename):
            os.remove(table_filename)

    if backend == "slurm":
        for job in jobs:
            _update_job_status(batch_filename, program_name, table, job["jobname"], status="submitting")
            job_id, stdout, stderr = _submit_slurm_command(
                job,
                cpus,
                time_estimate,
                memory,
                constraint=constraint,
            )
            status = "submitted" if job_id else "submit_failed"
            _update_job_status(
                batch_filename,
                program_name,
                table,
                job["jobname"],
                job_id=job_id,
                status=status,
                returncode="" if job_id else 1,
            )
            if stderr:
                text = stderr.decode(errors="replace") if isinstance(stderr, bytes) else str(stderr)
                if text.strip():
                    print(f"{job['jobname']} stderr: {text.strip()}")
        print(f"{len(jobs)} tasks added to slurm queue")
    else:
        for done, job in enumerate(jobs, start=1):
            _run_local_command(job, program_name, batch_filename, table)
            print_progress_bar(done, len(jobs), program_name=program_name, object_name="task(s)")

    _, table = _read_batch_table(batch_filename)
    if save_result_table:
        table_filename = _write_result_table(batch_filename, table)
        _save_result_table(table, table_filename)
    else:
        print(f"Logs folder: {os.path.dirname(jobs[0]['log_filename']) if jobs else os.path.dirname(batch_filename)}")
        print(f"Batch log: {batch_filename}")
    return table.sort_values(by="jobname").reset_index(drop=True)


def mixcr4_analyze_batch(sample_df, output_folder, command_template=None,
                         mixcr_path="mixcr", memory=32, time_estimate=1.5,
                         custom_tag_pattern_column=None, backend="local", cpus=40,
                         constraint=None):
    
    """
    Function for batch runs of MiXCR software.
    Runs commands locally by default and can also submit them to SLURM.
    By default this function uses `mixcr analyze` command for MiLab Hum RNA TCR Kit (with UMI). 
    To change the command template use `command_template` parameter

    Args:
        sample_df (pd.DataFrame): DataFrame, containing 'sample_id' column and 
            'R1' and 'R2' columns, containing paths (recommended full paths) to raw read files.
            Rows sharing the same R1/R2 pair are processed as a sample-barcoded mix and
            must also contain consistent 'mix_id' values and distinct ATGC-only 'SMPL' values.
            A non-empty 'miNNNPattern' value overrides the default sample-barcode tag pattern.
        output_folder (str): path to output folder
        command_template (str): MiXCR command template 
            (default: 'mixcr analyze milab-human-rna-tcr-umi-multiplex -f r1 r2 output_prefix').
            May be used as an example. Note that `mixcr analyze` and `r1 r2 output_prefix` are 
            "magical" parts of the template that should be kept as-is in the template, so change 
            only the part in-between these parts.
        mixcr_path (str): path to MiXCR binary
        memory (int): required OOM in GB
        time_estimate (numeric): time estimate in hours for the calculation. It
            is the limit for SLURM task
        backend (str): `local` or `slurm`
        cpus (int): CPU request for SLURM jobs
        constraint (str): Optional SLURM node constraint expression

    Returns:
        pd.DataFrame: submitted or completed job records
    """
    _validate_backend(backend)
    program_name="MIXCR4 Analyze Batch"

    # by default use the most popular preset for MiLaboratory Human TCR UMI MULTIPLEX Kit
    default_command_template = "mixcr analyze milab-human-rna-tcr-umi-multiplex -f r1 r2 output_prefix"
    if command_template is None:
        command_template = default_command_template

    required_columns = {"sample_id", "R1", "R2"}
    missing_columns = sorted(required_columns - set(sample_df.columns))
    if missing_columns:
        raise ValueError(f"sample_df is missing required columns: {', '.join(missing_columns)}")

    sample_df = sample_df.reset_index(drop=True)
    sample_barcoded_groups, sample_barcoded_positions = _sample_barcoded_groups(sample_df)
        
    # cut placeholders from command template
    command_template_parts = _strip_mixcr_template_placeholders(command_template)
    if sample_barcoded_groups and "--tag-pattern" in command_template_parts:
        raise ValueError("Please, remove '--tag-pattern' option from command_template for sample-barcoded files")

    # check input for custom tag pattern
    custom_tag_pattern = False
    if isinstance(custom_tag_pattern_column, str):
        if custom_tag_pattern_column not in sample_df.columns:
            raise ValueError(f"Specified tag-pattern columns '{custom_tag_pattern_column}' is not present in sample_df")
        if "--tag-pattern" in command_template_parts:
            raise ValueError(f"Please, remove '--tag-pattern' option from command_template, when you use custom tag-pattern")
        custom_tag_pattern = True
    
    # Create output dir if does not exist
    output_folder = os.path.abspath(output_folder)
    os.makedirs(output_folder, exist_ok=True)
    log_folder = os.path.join(output_folder, "logs")
    os.makedirs(log_folder, exist_ok=True)

    memory = _normalize_memory(memory)
        
    batch_filename = os.path.join(output_folder, "mixcr_analyze_batch.log")
    jobs = []
    
    # main cycle by samples and sample-barcoded mixes
    for i,r in sample_df.iterrows():
        if i in sample_barcoded_groups:
            sample_barcoded_group = sample_barcoded_groups[i]
            mix_id = sample_barcoded_group["mix_id"]
            barcode_table = sample_barcoded_group["rows"][["sample_id", "SMPL"]].rename(
                columns={"sample_id": "Sample"}
            )
            barcode_table["TagPattern"] = ""
            barcode_table = barcode_table[["Sample", "TagPattern", "SMPL"]]
            barcode_table_filename = os.path.join(output_folder, f"{mix_id}_samples.tsv")
            barcode_table.to_csv(barcode_table_filename, index=False, sep="\t")

            jobname = f"mixcr_analyze_{mix_id}"
            command = _mixcr_analyze_command(
                mixcr_path,
                memory,
                command_template_parts,
                r["R1"],
                r["R2"],
                mix_id,
                tag_pattern=sample_barcoded_group["tag_pattern"],
                extra_options=["--split-by-sample", "--sample-table", barcode_table_filename],
            )
            jobs.append({
                "jobname": jobname,
                "sample_id": mix_id,
                "command": command,
                "cwd": output_folder,
                "log_filename": os.path.join(log_folder, f"{jobname}.log"),
            })
            continue
        if i in sample_barcoded_positions:
            continue

        sample_id = r["sample_id"]
        r1 = r["R1"]
        r2 = r["R2"]
    #   output_prefix = os.path.join(output_folder, sample_id)
        output_prefix = sample_id
        tag_pattern = None
        if custom_tag_pattern:
            tag_pattern = r[custom_tag_pattern_column]
            if pd.isna(tag_pattern):
                raise ValueError(f"Empty tag pattern for sample '{sample_id}' in column '{custom_tag_pattern_column}'")
        command = _mixcr_analyze_command(
            mixcr_path,
            memory,
            command_template_parts,
            r1,
            r2,
            output_prefix,
            tag_pattern=tag_pattern,
        )
        jobname = f"mixcr_analyze_{sample_id}"
        jobs.append({
            "jobname": jobname,
            "sample_id": sample_id,
            "command": command,
            "cwd": output_folder,
            "log_filename": os.path.join(log_folder, f"{jobname}.log"),
        })

    return _run_mixcr_jobs(
        jobs,
        program_name,
        batch_filename,
        backend=backend,
        cpus=cpus,
        time_estimate=time_estimate,
        memory=memory,
        constraint=constraint,
    )


def _find_donor_clns_files(clns_filenames, sample_ids):
    sample_ids = {str(sample_id) for sample_id in sample_ids}
    matched_filenames = []
    for filename in clns_filenames:
        file_stem = os.path.basename(filename)[:-len(".clns")]
        if any(file_stem == sample_id or file_stem.endswith(f".{sample_id}")
               for sample_id in sample_ids):
            matched_filenames.append(filename)
    return matched_filenames


def find_alleles(input_dir, output_dir, mixcr_path="mixcr", sample_df=None,
                 backend="local", cpus=4, memory=32, time_estimate=0.5,
                 constraint=None):
    """
    Run MiXCR ``findAlleles`` once per donor and export the resulting clonotypes.

    Each donor job uses all ``.clns`` files in ``input_dir`` whose filename is
    either ``<sample_id>.clns`` or ends with ``.<sample_id>.clns``. The latter
    form supports clonosets produced by MiXCR sample splitting. After allele
    calling, every newly produced ``.clns`` file is exported to
    ``<filename>.clones.tsv`` in ``output_dir``.

    Args:
        input_dir (str): Folder containing preprocessed MiXCR ``.clns`` files.
        output_dir (str): Folder for allele-called ``.clns`` files, allele
            reports, libraries, exported clonotype tables, and job logs.
        mixcr_path (str): Path to the MiXCR binary.
        sample_df (pd.DataFrame): Sample metadata containing ``sample_id`` and
            ``donor_id`` columns. It may also be passed as the third positional
            argument when ``mixcr_path`` is omitted.
        backend (str): ``local`` or ``slurm``.
        cpus (int): CPU request for SLURM jobs.
        memory (int): MiXCR memory and SLURM memory request in GB.
        time_estimate (numeric): Time limit in hours for SLURM jobs.
        constraint (str): Optional SLURM node constraint expression.

    Returns:
        pd.DataFrame: Submitted or completed job records. Donors without any
        matching ``.clns`` files are printed and omitted from execution.
    """
    _validate_backend(backend)
    if sample_df is None and isinstance(mixcr_path, pd.DataFrame):
        sample_df = mixcr_path
        mixcr_path = "mixcr"
    if sample_df is None:
        raise ValueError("sample_df must be provided")
    required_columns = {"sample_id", "donor_id"}
    missing_columns = sorted(required_columns - set(sample_df.columns))
    if missing_columns:
        raise ValueError(f"sample_df is missing required columns: {', '.join(missing_columns)}")
    if sample_df[["sample_id", "donor_id"]].isna().any().any():
        raise ValueError("sample_df columns 'sample_id' and 'donor_id' must not contain empty values")
    if any(not str(value).strip()
           for value in sample_df[["sample_id", "donor_id"]].to_numpy().flat):
        raise ValueError("sample_df columns 'sample_id' and 'donor_id' must not contain empty values")

    input_dir = os.path.abspath(input_dir)
    if not os.path.isdir(input_dir):
        raise FileNotFoundError(f"Input directory does not exist: {input_dir}")
    output_dir = os.path.abspath(output_dir)
    os.makedirs(output_dir, exist_ok=True)
    log_folder = os.path.join(output_dir, "logs")
    os.makedirs(log_folder, exist_ok=True)

    memory = _normalize_memory(memory)
    clns_filenames = sorted(
        os.path.join(input_dir, filename)
        for filename in os.listdir(input_dir)
        if filename.endswith(".clns") and os.path.isfile(os.path.join(input_dir, filename))
    )

    jobs = []
    for donor_id, donor_samples in sample_df.groupby("donor_id", sort=False):
        donor_id = str(donor_id)
        sample_ids = donor_samples["sample_id"].astype(str).unique()
        donor_clns_filenames = _find_donor_clns_files(clns_filenames, sample_ids)

        print(f"Donor {donor_id}: {len(donor_clns_filenames)} .clns file(s)")
        for filename in donor_clns_filenames:
            print(f"  - {filename}")
        if not donor_clns_filenames:
            continue

        report_filename = os.path.join(output_dir, f"{donor_id}.findAlleles.report.txt")
        json_report_filename = os.path.join(output_dir, f"{donor_id}.findAlleles.report.json")
        mutations_filename = os.path.join(output_dir, f"{donor_id}_alleles.tsv")
        library_filename = os.path.join(output_dir, f"{donor_id}_alleles.json")
        output_template = os.path.join(output_dir, "{file_name}.clns")
        find_alleles_command = shlex.join([
            mixcr_path,
            f"-Xmx{memory}g",
            "findAlleles",
            "-f",
            "--report",
            report_filename,
            "--json-report",
            json_report_filename,
            "--export-alleles-mutations",
            mutations_filename,
            "--export-library",
            library_filename,
            "--output-template",
            output_template,
            *donor_clns_filenames,
        ])

        export_commands = []
        for input_filename in donor_clns_filenames:
            output_stem = os.path.basename(input_filename)[:-len(".clns")]
            output_clns_filename = os.path.join(output_dir, f"{output_stem}.clns")
            output_tsv_filename = os.path.join(output_dir, f"{output_stem}.clones.tsv")
            export_commands.append(shlex.join([
                mixcr_path,
                f"-Xmx{memory}g",
                "exportClones",
                output_clns_filename,
                output_tsv_filename,
            ]))

        jobname = f"mixcr_find_alleles_{donor_id}"
        jobs.append({
            "jobname": jobname,
            "sample_id": donor_id,
            "command": " && ".join([find_alleles_command, *export_commands]),
            "cwd": output_dir,
            "log_filename": os.path.join(log_folder, f"{jobname}.log"),
        })

    batch_filename = os.path.join(output_dir, "find_alleles_batch.log")
    return _run_mixcr_jobs(
        jobs,
        "MiXCR Find Alleles Batch",
        batch_filename,
        backend=backend,
        cpus=cpus,
        time_estimate=time_estimate,
        memory=memory,
        constraint=constraint,
    )


def find_shm_trees(input_dir, output_dir, mixcr_path="mixcr", sample_df=None,
                   backend="local", cpus=4, memory=32, time_estimate=0.5,
                   constraint=None):
    """
    Build and export MiXCR somatic hypermutation trees once per donor.

    Each donor job uses all allele-called ``.clns`` files in ``input_dir``
    whose filename is either ``<sample_id>.clns`` or ends with
    ``.<sample_id>.clns``. The job consecutively runs ``findShmTrees``,
    ``exportShmTreesWithNodes``, and ``exportShmTreesNewick``.

    Args:
        input_dir (str): Folder containing allele-called MiXCR ``.clns`` files.
        output_dir (str): Folder for ``.shmt``, TSV, Newick, report, and job log
            outputs.
        mixcr_path (str): Path to the MiXCR binary.
        sample_df (pd.DataFrame): Sample metadata containing ``sample_id`` and
            ``donor_id`` columns. It may also be passed as the third positional
            argument when ``mixcr_path`` is omitted.
        backend (str): ``local`` or ``slurm``.
        cpus (int): CPU request for SLURM jobs.
        memory (int): MiXCR memory and SLURM memory request in GB.
        time_estimate (numeric): Time limit in hours for SLURM jobs.
        constraint (str): Optional SLURM node constraint expression.

    Returns:
        pd.DataFrame: Submitted or completed job records. Donors without any
        matching ``.clns`` files are printed and omitted from execution.
    """
    _validate_backend(backend)
    if sample_df is None and isinstance(mixcr_path, pd.DataFrame):
        sample_df = mixcr_path
        mixcr_path = "mixcr"
    if sample_df is None:
        raise ValueError("sample_df must be provided")
    required_columns = {"sample_id", "donor_id"}
    missing_columns = sorted(required_columns - set(sample_df.columns))
    if missing_columns:
        raise ValueError(f"sample_df is missing required columns: {', '.join(missing_columns)}")
    if sample_df[["sample_id", "donor_id"]].isna().any().any():
        raise ValueError("sample_df columns 'sample_id' and 'donor_id' must not contain empty values")
    if any(not str(value).strip()
           for value in sample_df[["sample_id", "donor_id"]].to_numpy().flat):
        raise ValueError("sample_df columns 'sample_id' and 'donor_id' must not contain empty values")

    input_dir = os.path.abspath(input_dir)
    if not os.path.isdir(input_dir):
        raise FileNotFoundError(f"Input directory does not exist: {input_dir}")
    output_dir = os.path.abspath(output_dir)
    os.makedirs(output_dir, exist_ok=True)
    log_folder = os.path.join(output_dir, "logs")
    os.makedirs(log_folder, exist_ok=True)

    memory = _normalize_memory(memory)
    clns_filenames = sorted(
        os.path.join(input_dir, filename)
        for filename in os.listdir(input_dir)
        if filename.endswith(".clns") and os.path.isfile(os.path.join(input_dir, filename))
    )

    jobs = []
    for donor_id, donor_samples in sample_df.groupby("donor_id", sort=False):
        donor_id = str(donor_id)
        sample_ids = donor_samples["sample_id"].astype(str).unique()
        donor_clns_filenames = _find_donor_clns_files(clns_filenames, sample_ids)

        print(f"Donor {donor_id}: {len(donor_clns_filenames)} .clns file(s)")
        for filename in donor_clns_filenames:
            print(f"  - {filename}")
        if not donor_clns_filenames:
            continue

        trees_report = os.path.join(output_dir, f"{donor_id}_trees.log")
        output_shmt = os.path.join(output_dir, f"{donor_id}_trees.shmt")
        trees_export_filename = os.path.join(output_dir, f"{donor_id}_trees.tsv")
        trees_newick_dir = os.path.join(output_dir, f"{donor_id}_newick")
        find_trees_command = shlex.join([
            mixcr_path,
            f"-Xmx{memory}g",
            "findShmTrees",
            "-f",
            "--report",
            trees_report,
            *donor_clns_filenames,
            output_shmt,
        ])
        export_trees_command = shlex.join([
            mixcr_path,
            f"-Xmx{memory}g",
            "exportShmTreesWithNodes",
            output_shmt,
            trees_export_filename,
        ])
        export_newick_command = shlex.join([
            mixcr_path,
            f"-Xmx{memory}g",
            "exportShmTreesNewick",
            output_shmt,
            trees_newick_dir,
        ])

        jobname = f"mixcr_find_shm_trees_{donor_id}"
        jobs.append({
            "jobname": jobname,
            "sample_id": donor_id,
            "command": " && ".join([
                find_trees_command,
                export_trees_command,
                export_newick_command,
            ]),
            "cwd": output_dir,
            "log_filename": os.path.join(log_folder, f"{jobname}.log"),
        })

    batch_filename = os.path.join(output_dir, "find_shm_trees_batch.log")
    return _run_mixcr_jobs(
        jobs,
        "MiXCR Find SHM Trees Batch",
        batch_filename,
        backend=backend,
        cpus=cpus,
        time_estimate=time_estimate,
        memory=memory,
        constraint=constraint,
    )


def mixcr_7genes_run_batch(sample_df, output_folder, mixcr_path="mixcr", memory=32,
                           time_estimate=1.5, backend="local", cpus=40,
                           constraint=None):
    """
    Function for batch runs of the MiXCR software using the `mixcr analyze` command and the `Human 7GENES DNA Multiplex` MiXCR built-in preset.
    Incomplete rearrangements obtained by this kit are also included. For each incomplete rearrangement, unaligned reads from the previous 
    step are iteratively processed. Each output is stored in a subdirectory named after the corresponding rearrangement.
    Runs commands locally by default and can also submit them to SLURM.

    Args:
        sample_df (pd.DataFrame): DataFrame containing a 'sample_id' column and 
            'R1' and 'R2' columns containing paths (recommended full paths) to raw read files.
        output_folder (str): Path to the output folder.
        mixcr_path (str): Path to the MiXCR binary.
        memory (int): Required OOM in GB.
        time_estimate (numeric): Time estimate in hours for the calculation; it 
            is the limit for the SLURM task.
        backend (str): `local` or `slurm`.
        cpus (int): CPU request for SLURM jobs.
        constraint (str): Optional SLURM node constraint expression.

    Returns:
        pd.DataFrame: submitted or completed job records.
    """
    _validate_backend(backend)
    program_name="MIXCR4 Analyze 7genes Batch"
        
    # Create output dir if does not exist
    output_folder = os.path.abspath(output_folder)
    os.makedirs(output_folder, exist_ok=True)
    log_folder = os.path.join(output_folder, "logs")
    os.makedirs(log_folder, exist_ok=True)

    memory = _normalize_memory(memory)
        
    batch_filename = os.path.join(output_folder, "mixcr_analyze_7genes_batch.log")
    
    list_of_incomplete_rearrangements = ["DJ_TRB", "VDD_TRD", "DDJ_TRD", "DD_TRD", "DJ_IGH", "VKDE_IGK", "CINTRON_KDE_IGK"]

    jobs = []
    # main cycle by samples
    for i,r in sample_df.iterrows():
        sample_id = r["sample_id"]
        r1 = r["R1"]
        r2 = r["R2"]
        output_prefix = sample_id
        
        R1na = f"{sample_id}_R1_not_aligned.fastq.gz"
        R2na = f"{sample_id}_R2_not_aligned.fastq.gz"
        
        commands = []
        first_command = f'{mixcr_path} -Xmx{memory}g analyze milab-human-dna-xcr-7genes-multiplex -f --not-aligned-R1 {R1na} --not-aligned-R2 {R2na} {r1} {r2} {output_prefix}'
        commands.append(first_command)
        
        for rearrangement in list_of_incomplete_rearrangements:
            
            # swap r and Rna so we would not implement copy of R_na
            r1, R1na = R1na, r1
            r2, R2na = R2na, r2
            
            output_prefix = os.path.join(rearrangement, sample_id)
            
            R1na = f"{output_prefix}_R1_not_aligned.fastq.gz"
            R2na = f"{output_prefix}_R2_not_aligned.fastq.gz"
            
            i_r_command = f'{mixcr_path} -Xmx{memory}g analyze generic-amplicon -f --species hs --library {rearrangement} --assemble-clonotypes-by CDR3 --dna --floating-left-alignment-boundary --floating-right-alignment-boundary J -MexportClones.splitFilesBy=[] --not-aligned-R1 {R1na} --not-aligned-R2 {R2na} {r1} {r2} {output_prefix}'
            commands.append(i_r_command)
        
        jobname = f"mixcr_analyze_{sample_id}"
        
        # join commands by && so that next command runs if previous was finished without error and add new lines to the script
        command = " && \\ \n".join(commands)
        jobs.append({
            "jobname": jobname,
            "sample_id": sample_id,
            "command": command,
            "cwd": output_folder,
            "log_filename": os.path.join(log_folder, f"{jobname}.log"),
        })

    return _run_mixcr_jobs(
        jobs,
        program_name,
        batch_filename,
        backend=backend,
        cpus=cpus,
        time_estimate=time_estimate,
        memory=memory,
        constraint=constraint,
    )


def mixcr4_reports(folder, mixcr_path="mixcr", backend="local",
                   cpus=40, time_estimate=1, memory=32, constraint=None):
    
    """
    runs `mixcr exportQc` commands - `align`, `chainUsage` and `tags` in a given folder 
    for all `.clns` filenames. `align` and `chainUsage` are run twice to create both 
    `svg` and `pdf` files.

    Args:
        folder (str): folder in which to run the `mixcr exportQc` commands
        mixcr_path (str): path to MiXCR binary
        backend (str): `local` or `slurm`
        cpus (int): CPU request for SLURM jobs
        time_estimate (numeric): time estimate in hours for SLURM jobs
        memory (int): MiXCR memory in GB
        constraint (str): Optional SLURM node constraint expression
    Returns:
        pd.DataFrame: submitted or completed job records

    """
    _validate_backend(backend)

    program_name="MiXCR 4 Reports"
    memory = _normalize_memory(memory)
    folder = os.path.abspath(folder)
    os.makedirs(folder, exist_ok=True)
    log_folder = os.path.join(folder, "logs")
    os.makedirs(log_folder, exist_ok=True)
    
    # clns_filenames = os.path.join(folder, "*.clns")
    # align_filename = os.path.join(folder, "alignQc.png")
    # chains_filename = os.path.join(folder, "chainsQc.png")
    # tags_filename = os.path.join(folder, "tagsQc.pdf")
    clns_filenames = "*.clns"
    align_filename = "alignQc.svg"
    chains_filename = "chainsQc.svg"
    align_filename_pdf = "alignQc.pdf"
    chains_filename_pdf = "chainsQc.pdf"
    tags_filename = "tagsQc.pdf"
    #tables_filename = os.path.join(folder, "tables.tsv")
    #preproc_filename = os.path.join(folder, "preproc_tables.tsv")
    #postanalysis_filename = os.path.join(folder, "postanalysis.json")
    

    
    commands = {"alignQc": f"{mixcr_path} -Xmx{memory}g exportQc align -f {clns_filenames} {align_filename}",
                "chainUsage": f"{mixcr_path} -Xmx{memory}g exportQc chainUsage -f {clns_filenames} {chains_filename}",
                "alignQcPDF": f"{mixcr_path} -Xmx{memory}g exportQc align -f {clns_filenames} {align_filename_pdf}",
                "chainUsagePDF": f"{mixcr_path} -Xmx{memory}g exportQc chainUsage -f {clns_filenames} {chains_filename_pdf}",
                "tagsQc": f"{mixcr_path} -Xmx{memory}g exportQc tags -f {clns_filenames} {tags_filename}"#,
                #"postanalysis": f"{mixcr_path} -Xmx32g postanalysis individual -f --default-downsampling none --default-weight-function umi --only-productive --tables {tables_filename} --preproc-tables {preproc_filename} {clns_filenames} {postanalysis_filename}"
               }
    

    batch_filename = os.path.join(folder, "mixcr_reports_slurm_batch.log")
    jobs = []
    for jobname, command in commands.items():
        jobs.append({
            "jobname": jobname,
            "sample_id": "",
            "command": command,
            "cwd": folder,
            "log_filename": os.path.join(log_folder, f"{jobname}.log"),
        })

    return _run_mixcr_jobs(
        jobs,
        program_name,
        batch_filename,
        backend=backend,
        cpus=cpus,
        time_estimate=time_estimate,
        memory=memory,
        save_result_table=False,
        constraint=constraint,
    )


def get_processing_table(folder, show_offtarget=False, offtarget_chain_threshold=0.01,
                         small=False):
    """
    Searches for clonosets in the the folder, extracts their sample_id's and shows main
    processing stats in a table format. By default does not show "off-target" clonosets - 
    those having less than 1% (default, may be overriden) of reads for the sample_id.
    For example, you have sequenced TRB sample, but there is found 0.5% (by read count) 
    of TRA chains for the same sample_id, then the clonoset will not be shown in the table.
    You can specify `show_offtarget=True` to display all found chains in the table or 
    outherwise set a higher value for `offtarget_chain_threshold` (`0.01` by default).

    Args:
        folder (str or list): folder or list of folders in which to look for clonosets and
            processing stats
        show_offtarget (bool): add offtarget chains to the stats
        offtarget_chain_threshold (float): threshold for off-target chains
        small (bool): return only the main processing stats columns
    
    Returns:
        df (pd.DataFrame): dataframe, containing `sample_id`, `extracted_chain` and
            different processing stats columns. If clonoset filenames contain both
            `mix_id` and `sample_id`, the dataframe also contains a `mix_id` column.
            There may be several rows with the same `sample_id`, with each found
            `extracted_chain`.
    """

    if isinstance(folder, list):
        tables = []
        for f in folder:
            table = get_processing_table(
                f,
                show_offtarget=show_offtarget,
                offtarget_chain_threshold=offtarget_chain_threshold,
                small=small,
            )
            tables.append(table)
        result_df = pd.concat(tables)
        if "mix_id" in result_df.columns:
            columns = ["sample_id", "mix_id"] + [
                column for column in result_df.columns
                if column not in {"sample_id", "mix_id"}
            ]
            result_df = result_df[columns]
        return result_df.sort_values(by="sample_id").reset_index(drop=True)
    
    results = []
    clonosets = find_all_exported_clonosets_in_folder(folder, chain=None)
    has_mix_id = "mix_id" in clonosets.columns

    for i, r in clonosets.iterrows():
        sample_id = r["sample_id"]
        mix_id = r["mix_id"] if has_mix_id else None
        report_sample_id = sample_id
        align_report_id = sample_id
        if pd.notna(mix_id):
            report_sample_id = f"{mix_id}.{sample_id}"
            align_report_id = mix_id
        chain = r["chain"]
        align_report = read_json_report(align_report_id, folder, "align")
        
        try:
            refine_report = read_json_report(report_sample_id, folder, "refine")
            umi = True
        except FileNotFoundError:
            umi = False
            
        assemble_report = read_json_report(report_sample_id, folder, "assemble")

        # print(sample_id, chain)
        clonoset = read_clonoset(r.filename)
        clonoset_f = filter_by_functionality(clonoset)

        # align report
        Rt=align_report["totalReadsProcessed"]
        Ru=align_report["totalReadsProcessed"]-align_report["notAlignedReasons"]["NoBarcode"]
        Ru_pc = round(Ru/Rt*100, 2)
        Ra=align_report["aligned"]
        Ra_pc = round(Ra/Rt*100, 2)
        Roa = align_report["overlappedAligned"]
        Roa_pc = round(Roa/Ra*100, 2)
        
        if umi:
        #Ra2=refine_report["correctionReport"]["inputRecords"] ##### differs from Ra, but D.Bolotin did not explain why
        
            UMIa=refine_report["correctionReport"]["steps"][0]["inputDiversity"]
            UMIc=refine_report["correctionReport"]["steps"][0]["outputDiversity"]
            try:
                UMIf=refine_report["correctionReport"]["filterReport"]["numberOfGroupsAccepted"]
            except TypeError:
                UMIf=UMIc
            Rf=refine_report["correctionReport"]["outputRecords"]
            try:
                overseq_threshold = int(refine_report["correctionReport"]["filterReport"]["operatorReports"][0]["operatorReport"]["threshold"])
            except TypeError:
                overseq_threshold = None
            reads_per_umi = round(Rf/UMIf, 2)
        else:
            UMIa = np.nan
            UMIc = np.nan
            UMIf = np.nan
            Rf = np.nan
            overseq_threshold = np.nan
            reads_per_umi = np.nan
            
        Ct=assemble_report["clones"]
        Rcl=assemble_report["readsInClones"]
        
        Ctc=len(clonoset)
        Rclc=int(clonoset.readCount.sum())
        
        Cfunc=len(clonoset_f)
        Rfunc=int(clonoset_f.readCount.sum())
        if umi:
            UMIcl=clonoset.uniqueMoleculeCount.sum()
            UMIfunc=clonoset_f.uniqueMoleculeCount.sum()
        else:
            UMIcl=np.nan
            UMIfunc=np.nan
        if umi and overseq_threshold is None:
            reads_per_umi = round(Rclc/UMIcl, 2)

        result = [sample_id, chain, Rt, Ru_pc, Ra_pc, Roa_pc, UMIa, UMIc, overseq_threshold, Rf, UMIf, reads_per_umi, Ct, Rcl, Ctc, Rclc, Cfunc, Rfunc, UMIcl, UMIfunc]
        if has_mix_id:
            result.insert(1, mix_id)
        results.append(result)
    result_columns = ["sample_id", "extracted_chain", "reads_total", "reads_with_umi_pc", "reads_aligned_pc", "reads_overlapped_aln_pc",
                      "total_umi", "umi_after_correction", "overseq_threshold", "reads_after_filter", "umi_after_filter",
                      "reads_per_umi", "clones_total", "reads_in_clones_total", "clones", "reads_in_clones", "clones_func", "reads_in_func_clones", "umi_in_clones", "umi_in_func_clones"]
    if has_mix_id:
        result_columns.insert(1, "mix_id")
    result_df = pd.DataFrame(results, columns=result_columns)
    if not show_offtarget:
        result_df = result_df.loc[result_df.reads_in_clones/result_df.reads_in_clones_total > offtarget_chain_threshold]
    if small:
        small_columns = [
            "sample_id", "extracted_chain", "reads_total", "reads_with_umi_pc",
            "reads_aligned_pc", "reads_overlapped_aln_pc", "reads_per_umi",
            "overseq_threshold", "clones_func", "umi_in_func_clones",
        ]
        if has_mix_id:
            small_columns.insert(1, "mix_id")
        result_df = result_df[small_columns]
    return result_df.sort_values(by="sample_id").reset_index(drop=True)


def show_report_images(folder):
    """
    Display MiXCR QC report images in a Jupyter notebook.

    The function looks for `alignQc.svg` and `chainsQc.svg` in `folder`.
    If an SVG file is missing, it falls back to the corresponding PNG file:
    `alignQc.png` or `chainsQc.png`. If neither image exists for a report,
    a short message is printed and execution continues.

    These images can be generated with `mixcr4_reports`.

    Args:
        folder (str): folder in which to look for MiXCR JSON reports.
    
    Returns:
        None.

    """
    try:
        from IPython.display import Image, display, SVG
    except ImportError as exc:
        raise ImportError(
            "Displaying MiXCR report images requires IPython. "
            "Install it with `pip install ipython`."
        ) from exc
    
    svg_align_filename = os.path.join(folder, "alignQc.svg")
    svg_chain_filename = os.path.join(folder, "chainsQc.svg")
    png_align_filename = os.path.join(folder, "alignQc.png")
    png_chain_filename = os.path.join(folder, "chainsQc.png")
    
    if os.path.exists(svg_align_filename):
        print(svg_align_filename)
        display(SVG(filename=svg_align_filename))
    elif os.path.exists(png_align_filename):
        print(png_align_filename)
        display(Image(filename=png_align_filename))
    else:
        print("No alignQc image found (svg or png)")

    if os.path.exists(svg_chain_filename):
        print(svg_chain_filename)
        display(SVG(filename=svg_chain_filename))
    elif os.path.exists(png_chain_filename):
        print(png_chain_filename)
        display(Image(filename=png_chain_filename))
    else:
        print("No chainQc image found (svg or png)")


def _qc_plot_sample_labels(table):
    labels = table["sample_id"].astype(str)
    if "extracted_chain" in table.columns and labels.duplicated(keep=False).any():
        labels = labels + " (" + table["extracted_chain"].astype(str) + ")"
    return labels


def _coverage_table_from_refine_reports(folder):
    rows = []
    try:
        filenames = os.listdir(folder)
    except FileNotFoundError:
        print("No such file or directory")
        return pd.DataFrame(columns=["sample_id", "reads_per_umi", "overseq_threshold"])

    for filename in filenames:
        match = re.match(r"(.+)\.refine\.report\.json$", filename)
        if match is None:
            continue
        sample_id = match.group(1)
        report = open_json_report(os.path.join(folder, filename))
        correction_report = report.get("correctionReport", {})
        reads_after_filter = correction_report.get("outputRecords")
        umi_after_filter = None
        overseq_threshold = np.nan
        filter_report = correction_report.get("filterReport")
        if isinstance(filter_report, dict):
            umi_after_filter = filter_report.get("numberOfGroupsAccepted")
            try:
                overseq_threshold = int(
                    filter_report["operatorReports"][0]["operatorReport"]["threshold"]
                )
            except (KeyError, IndexError, TypeError, ValueError):
                overseq_threshold = np.nan
        if umi_after_filter is None:
            try:
                umi_after_filter = correction_report["steps"][0]["outputDiversity"]
            except (KeyError, IndexError, TypeError):
                umi_after_filter = None
        if reads_after_filter is None or umi_after_filter in [None, 0]:
            reads_per_umi = np.nan
        else:
            reads_per_umi = round(reads_after_filter / umi_after_filter, 2)
        rows.append({
            "sample_id": sample_id,
            "reads_per_umi": reads_per_umi,
            "overseq_threshold": overseq_threshold,
        })
    return pd.DataFrame(rows, columns=["sample_id", "reads_per_umi", "overseq_threshold"])


def _format_qc_value(value):
    numeric_value = float(value)
    if numeric_value.is_integer():
        return str(int(numeric_value))
    return str(value)


def _plot_coverage_qc(folder, processing_table=None, output_file=None,
                      show_offtarget=False, offtarget_chain_threshold=0.01):
    del show_offtarget, offtarget_chain_threshold
    if processing_table is None:
        processing_table = _coverage_table_from_refine_reports(folder)
    required_columns = {"sample_id", "reads_per_umi", "overseq_threshold"}
    missing = required_columns - set(processing_table.columns)
    if missing:
        missing_text = ", ".join(sorted(missing))
        raise ValueError(f"Coverage plot requires columns: {missing_text}")

    plot_data = processing_table.copy()
    plot_data = plot_data.dropna(subset=["reads_per_umi"])
    if len(plot_data) == 0:
        print("No coverage data found")
        return
    plot_data = plot_data.sort_values(by="sample_id", ascending=False).reset_index(drop=True)
    plot_data.index = _qc_plot_sample_labels(plot_data)

    size = plot_data.shape[0]
    bar_height = 0.85
    min_size = 7
    min_size_2 = 10
    plot_rows = max(size, min_size)
    if size > min_size:
        plot_rows = max(size, min_size_2)
    fig, ax = plt.subplots(figsize=(9, plot_rows * bar_height * 0.5), dpi=100, constrained_layout=True)
    y = np.arange(len(plot_data))
    bars = ax.barh(
        y=y,
        width=plot_data["reads_per_umi"].values,
        height=bar_height,
        color="#d8c3a5",
        label="Reads per UMI",
    )
    for bar, reads_per_umi in zip(bars, plot_data["reads_per_umi"]):
        ax.annotate(
            _format_qc_value(reads_per_umi),
            xy=(bar.get_width(), bar.get_y() + bar.get_height() / 2),
            xytext=(4, 0),
            textcoords="offset points",
            ha="left",
            va="center",
            color="#5c5142",
        )
    marker_labeled = False
    for i, threshold in enumerate(plot_data["overseq_threshold"]):
        if pd.isna(threshold):
            continue
        marker_position = float(threshold) - 1
        ax.vlines(
            marker_position,
            i - bar_height / 2,
            i + bar_height / 2,
            color="#d62728",
            linewidth=2,
            label="Overseq threshold - 1" if not marker_labeled else None,
        )
        ax.annotate(
            _format_qc_value(threshold),
            xy=(marker_position, i),
            xytext=(4, 0),
            textcoords="offset points",
            ha="left",
            va="center",
            color="#d62728",
        )
        marker_labeled = True
    marker_positions = plot_data["overseq_threshold"].dropna().astype(float) - 1
    rightmost_value = max(
        float(plot_data["reads_per_umi"].max()),
        float(marker_positions.max()) if not marker_positions.empty else 0,
    )
    ax.set_xlim(0, rightmost_value + max(1, rightmost_value * 0.12))
    ax.set_yticks(y)
    ax.set_yticklabels(plot_data.index)
    ax.set_ylim(-0.5, len(plot_data) - 0.5)
    ax.set_xlabel("reads per UMI")
    fig.legend(loc="outside upper center", title="Coverage", ncol=2, frameon=False)
    sns.despine(left=True, bottom=True)
    plt.show()
    if output_file is not None:
        ax.get_figure().savefig(output_file, bbox_inches="tight")


def show_qc_plot(folder, chart_type='align', count_type='percent', output_file=None,
                 processing_table=None, show_offtarget=False,
                 offtarget_chain_threshold=0.01):
    """
    Shows quality control reports in MiXCR-like style

    Args:
        folder (str): folder in which to look for QC images.
        chart_type (str): Possible values are `align` (corresponds to
            `mixcr exportQc align`), `chains` (plots `clonalChainUsage` from
            `*.assemble.report.json` files), or `coverage` (plots
            `reads_per_umi` and `overseq_threshold - 1` directly from
            `*.refine.report.json` files). Deprecated alias `summary` is
            accepted as `align`.
        count_type (str): possible values are: `percent`, `abs`
        output_file (str): filename ending with '.png' to save an output plot to
        processing_table (pd.DataFrame): optional precomputed processing table
            for `chart_type="coverage"`.
        show_offtarget (bool): ignored for coverage plots, kept for backwards
            compatibility.
        offtarget_chain_threshold (float): ignored for coverage plots, kept for
            backwards compatibility.

    Returns:
        None

    """
    if chart_type == "summary":
        warnings.warn(
            "`chart_type='summary'` is deprecated; use `chart_type='align'` instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        chart_type = "align"
    if chart_type not in ["align", "chains", "coverage"]:
        raise ValueError("chart_type must be one of: 'align', 'chains', 'coverage'")
    if count_type not in ["percent", "abs"]:
        raise ValueError("count_type must be either 'percent' or 'abs'")
    if chart_type == "coverage":
        return _plot_coverage_qc(
            folder,
            processing_table=processing_table,
            output_file=output_file,
            show_offtarget=show_offtarget,
            offtarget_chain_threshold=offtarget_chain_threshold,
        )

    CHAIN_VARIANTS = ['IGH', 'IGK', 'IGL', 'TRA', 'TRB', 'TRD', 'TRG']
    files = []
    try:
        all_files = os.listdir(folder)
        files = []
        for f in all_files:
            match = re.match(r"(.+)\.([a-zA-Z0-9_]+)\.report\.json", f)
            if match is not None:
                sample_id = match.group(1)
                report_type = match.group(2)
                files.append([sample_id, report_type])
            else:
                continue
    except FileNotFoundError:
        print('No such file or directory')
    df_list = []
    expected_report_type = 'align' if chart_type == 'align' else 'assemble'
    for file in files:
        report_type = file[1]
        if report_type == expected_report_type:
            json_report_contents = read_json_report(file[0], folder, report_type=report_type)
            
            if chart_type == 'align':
                align_data = json_report_contents['notAlignedReasons']
                renaming_dict = {'NoHits': 'No hits (not TCR/IG?)',
                                'NoCDR3Parts': 'No CDR3 parts',
                                'NoVHits': 'No V hits',
                                'NoJHits': 'No J hits',
                                'VAndJOnDifferentTargets': 'No target with both V and J',
                                'LowTotalScore': 'Low total score',
                                'NoBarcode': 'Absent barcode',
                                'SampleNotMatched': 'Sample not matched',
                                }
                align_df = {}
                for old, new in renaming_dict.items():
                    align_df[new] = [align_data.get(old, 0)]
                align_df['Successfully aligned'] = json_report_contents['aligned']
                df_list.append(pd.DataFrame(align_df, index=[file[0]])) 
                
            elif chart_type == 'chains':
                chain_usage_data = json_report_contents['clonalChainUsage']['chains']
                align_df = {}
                for chain, data in chain_usage_data.items(): 
                    align_df.update({chain: data['total'] - data['nonFunctional'],
                                f'{chain} (OOF)': data['isOOF'],
                                f'{chain} (stops)': data['hasStops']})
                df_list.append(pd.DataFrame(align_df, index=[file[0]]))
    if len(df_list) == 0:
        print(f"No {chart_type} report data found")
        return
    results = pd.concat(df_list)
    results = results.sort_index(ascending=False)
    results = results.fillna(0)
    if count_type == 'percent':
        results =  results.div(results.sum(axis=1), axis=0) * 100
    if chart_type == 'align':
        order = ['Successfully aligned', 
                 'No hits (not TCR/IG?)', 
                 'No CDR3 parts', 
                 'No V hits', 
                 'No J hits', 
                 'No target with both V and J', 
                 'Low total score', 
                 'Absent barcode'] 
        colors = ['#3ecd8d', '#fed470', '#fda163', '#f36c5a', '#d64470', '#a03080', '#702084', '#451777']
        colormap = ListedColormap(colors=colors,
                                    name='mixcr')
    elif chart_type == 'chains':
        column_chains = sorted({column.split(" ")[0] for column in results.columns})
        chains_found = [chain for chain in CHAIN_VARIANTS if chain in column_chains]
        chains_found += [chain for chain in column_chains if chain not in chains_found]
        order = []
        for chain in chains_found:
            order += [
                column for column in [chain, f"{chain} (stops)", f"{chain} (OOF)"]
                if column in results.columns
            ]
        chain_color_map = {
            'IGH': ('#c26a27', '#ffcb8f', '#ff9429'),
            'IGK': ('#a324b2', '#faaafa', '#e553e5'),
            'IGL': ('#ad3757', '#ffadba', '#f05670'),
            'TRA': ('#105bcc', '#99ccff', '#2d93fa'),
            'TRB': ('#198020', '#99e099', '#42b842'),
            'TRD': ('#068a94', '#90e0e0', '#27c2c2'),
            'TRG': ('#5f31cc', '#c1adff', '#845cff'),
        }
        colors = []
        for column in order:
            chain = column.split(" ")[0]
            color_index = 0
            if "(OOF)" in column:
                color_index = 1
            elif "(stops)" in column:
                color_index = 2
            colors.append(chain_color_map.get(chain, ('#808080', '#c0c0c0', '#a0a0a0'))[color_index])
    results = results[order]
    size = results.shape[0]
    # ax = results.plot.barh(width=0.85, figsize=(9, size * 0.5),  stacked=True, colormap=colormap)
    bar_height = 0.85
    min_size = 7
    min_size_2 = 10
    plot_rows = max(size, min_size)
    if size > min_size:
            plot_rows = max(size, min_size_2)
    fig, ax = plt.subplots(figsize=(9, plot_rows * bar_height * 0.5), dpi=100, constrained_layout=True)
    ax.set_ylim(-0.5, results.shape[0] - 0.5)
    bottom = np.zeros(len(results))
    for i, column in enumerate(results.columns):
        values = results[column].values
        ax.barh(
            y=results.index,
            width=values,
            height=bar_height,
            left=bottom,
            label=column,
            color=colors[i],
        )
        bottom = [b + v for b, v in zip(bottom, values)]
    if count_type == 'percent':
        ax.set_xlabel('%')
    else:
        ax.set_xlabel('read count')
    if chart_type == 'align':
        fig.legend(loc='outside upper center',  title='Alignments rate', ncol=2, frameon=False)
    elif chart_type == 'chains':
        fig.legend(loc='outside upper center',  title='Clonal chain usage', ncol=max(1, len(chains_found)), frameon=False)
    # plt.tight_layout()
    sns.despine(left=True, bottom=True)
    plt.show()
    if output_file is not None:
        ax.get_figure().savefig(output_file, bbox_inches='tight')
    return


def show_report_images_new(*args, **kwargs):
    """
    Deprecated alias for `show_qc_plot`.

    Use `show_qc_plot` instead. This alias will be removed in a future version.
    """
    warnings.warn(
        "`show_report_images_new` is deprecated; use `show_qc_plot` instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    return show_qc_plot(*args, **kwargs)
