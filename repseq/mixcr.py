from .common_functions import print_progress_bar, filter_by_functionality
import pandas as pd
import numpy as np
import concurrent.futures
import os
import shlex
import subprocess
from time import sleep
from .slurm import run_slurm_command_from_jupyter
from .io import read_json_report, read_clonoset
from .clonosets import find_all_exported_clonosets_in_folder
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import seaborn as sns
import re


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


def create_batch_file(filename, program_name, tasks_num, backend="local"):
    with open(filename, "w") as f:
        f.write(f"# {program_name}: {backend} run {tasks_num} tasks\n")


def _append_batch_status(filename, jobname, status):
    with open(filename, "a") as f:
        f.write(f"{jobname} {status}\n")


def check_batch_progress(filename, loop=False):
    with open(filename, "r") as f:
        first_line = f.readlines()[0]
    match = re.match(r'# ([^:]+): ([^ ]+) run ([0-9]+) tasks', first_line)
    if match is None:
        raise ValueError("Could not parse batch progress file header")
    program_name = match[1]
    tasks = int(match[3])

    while True:
        with open(filename, "r") as f:
            text = f.read()
        finished = len(re.findall(r"\bfinished\b", text))
        failed = len(re.findall(r"\bfailed\b", text))
        if not loop:
            print(text)
        print_progress_bar(finished + failed, tasks, program_name=program_name, object_name="task(s)")
        if not loop or finished + failed == tasks:
            break
        sleep(0.5)


def _run_local_command(command, jobname, cwd, log_filename, batch_filename):
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
    _append_batch_status(batch_filename, jobname, status)
    return {
        "jobname": jobname,
        "backend": "local",
        "command": command,
        "cwd": cwd,
        "log_filename": log_filename,
        "returncode": process.returncode,
        "status": status,
        "stdout": "",
        "stderr": "",
    }


def _submit_slurm_command(command, jobname, cwd, cpus, time_estimate, memory, log_filename, batch_filename):
    slurm_command = command
    if cwd is not None:
        slurm_command = f"cd {shlex.quote(cwd)} && {command}"
    slurm_command += f'; echo "{jobname} finished" >> {shlex.quote(batch_filename)}'
    stdout, stderr = run_slurm_command_from_jupyter(
        slurm_command,
        jobname,
        cpus,
        time_estimate,
        memory,
        log_filename=log_filename,
    )
    return {
        "jobname": jobname,
        "backend": "slurm",
        "command": slurm_command,
        "cwd": cwd,
        "log_filename": log_filename,
        "returncode": None,
        "status": "submitted",
        "stdout": stdout.decode(errors="replace") if isinstance(stdout, bytes) else stdout,
        "stderr": stderr.decode(errors="replace") if isinstance(stderr, bytes) else stderr,
    }


def _run_mixcr_jobs(jobs, program_name, batch_filename, backend="local", max_workers=1,
                    cpus=40, time_estimate=1.5, memory=32):
    _validate_backend(backend)
    create_batch_file(batch_filename, program_name, len(jobs), backend=backend)

    if backend == "slurm":
        results = []
        for job in jobs:
            result = _submit_slurm_command(
                job["command"],
                job["jobname"],
                job["cwd"],
                cpus,
                time_estimate,
                memory,
                job["log_filename"],
                batch_filename,
            )
            print(result["jobname"], result["stdout"], result["stderr"])
            results.append(result)
        print(f"{len(jobs)} tasks added to slurm queue\n")
    else:
        if not isinstance(max_workers, int) or max_workers < 1:
            raise ValueError("max_workers must be a positive integer")
        print(f"Running {len(jobs)} tasks locally with max_workers={max_workers}\n")
        if max_workers == 1:
            results = [
                _run_local_command(job["command"], job["jobname"], job["cwd"], job["log_filename"], batch_filename)
                for job in jobs
            ]
        else:
            with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
                futures = [
                    executor.submit(
                        _run_local_command,
                        job["command"],
                        job["jobname"],
                        job["cwd"],
                        job["log_filename"],
                        batch_filename,
                    )
                    for job in jobs
                ]
                results = [future.result() for future in concurrent.futures.as_completed(futures)]

    print(f'To see progress, run:\nmixcr.check_batch_progress("{batch_filename}", loop=True)')
    print(f'To see current progress:\nmixcr.check_batch_progress("{batch_filename}")')
    columns = ["jobname", "backend", "command", "cwd", "log_filename", "returncode", "status", "stdout", "stderr"]
    if len(results) == 0:
        return pd.DataFrame(columns=columns)
    return pd.DataFrame(results, columns=columns).sort_values(by="jobname").reset_index(drop=True)


def mixcr4_analyze_batch(sample_df, output_folder, command_template=None,
                         mixcr_path="mixcr", memory=32, time_estimate=1.5,
                         custom_tag_pattern_column=None, backend="local",
                         max_workers=1, cpus=40):
    
    """
    Function for batch runs of MiXCR software.
    Runs commands locally by default and can also submit them to SLURM.
    By default this function uses `mixcr analyze` command for MiLab Hum RNA TCR Kit (with UMI). 
    To change the command template use `command_template` parameter

    Args:
        sample_df (pd.DataFrame): DataFrame, containing 'sample_id' column and 
            'R1' and 'R2' columns, containing paths (recommended full paths) to raw read files
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
        max_workers (int): number of local sample jobs to run in parallel
        cpus (int): CPU request for SLURM jobs

    Returns:
        pd.DataFrame: submitted or completed job records
    """
    _validate_backend(backend)
    program_name="MIXCR4 Analyze Batch"

    # by default use the most popular preset for MiLaboratory Human TCR UMI MULTIPLEX Kit
    default_command_template = "mixcr analyze milab-human-rna-tcr-umi-multiplex -f r1 r2 output_prefix"
    if command_template is None:
        command_template = default_command_template
        
    # cut placeholders from command template
    remove_list = ["mixcr", "r1", "r2", "output_prefix"]
    command_template = ' '.join([w for w in command_template.split() if w not in remove_list])

    # check input for custom tag pattern
    custom_tag_pattern = False
    if isinstance(custom_tag_pattern_column, str):
        if custom_tag_pattern_column not in sample_df.columns:
            raise ValueError(f"Specified tag-pattern columns '{custom_tag_pattern_column}' is not present in sample_df")
        if "--tag-pattern" in command_template.split():
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
    
    # main cycle by samples
    for i,r in sample_df.iterrows():
        sample_id = r["sample_id"]
        r1 = r["R1"]
        r2 = r["R2"]
    #   output_prefix = os.path.join(output_folder, sample_id)
        output_prefix = sample_id
        if custom_tag_pattern:
            tag_pattern = r[custom_tag_pattern_column]
            command = f'{mixcr_path} -Xmx{memory}g {command_template} --tag-pattern "{tag_pattern}" {r1} {r2} {output_prefix}'
        else:
            command = f'{mixcr_path} -Xmx{memory}g {command_template} {r1} {r2} {output_prefix}'
        jobname = f"mixcr_analyze_{sample_id}"
        jobs.append({
            "jobname": jobname,
            "command": command,
            "cwd": output_folder,
            "log_filename": os.path.join(log_folder, f"{jobname}.log"),
        })

    return _run_mixcr_jobs(
        jobs,
        program_name,
        batch_filename,
        backend=backend,
        max_workers=max_workers,
        cpus=cpus,
        time_estimate=time_estimate,
        memory=memory,
    )


def mixcr_7genes_run_batch(sample_df, output_folder, mixcr_path="mixcr", memory=32,
                           time_estimate=1.5, backend="local", max_workers=1, cpus=40):
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
        max_workers (int): number of local sample jobs to run in parallel.
        cpus (int): CPU request for SLURM jobs.

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
            "command": command,
            "cwd": output_folder,
            "log_filename": os.path.join(log_folder, f"{jobname}.log"),
        })

    return _run_mixcr_jobs(
        jobs,
        program_name,
        batch_filename,
        backend=backend,
        max_workers=max_workers,
        cpus=cpus,
        time_estimate=time_estimate,
        memory=memory,
    )


def mixcr4_reports(folder, mixcr_path="mixcr", backend="local", max_workers=1,
                   cpus=40, time_estimate=1, memory=32):
    
    """
    runs `mixcr exportQc` commands - `align`, `chainUsage` and `tags` in a given folder 
    for all `.clns` filenames. `align` and `chainUsage` are run twice to create both 
    `svg` and `pdf` files.

    Args:
        folder (str): folder in which to run the `mixcr exportQc` commands
        mixcr_path (str): path to MiXCR binary
        backend (str): `local` or `slurm`
        max_workers (int): number of local report jobs to run in parallel
        cpus (int): CPU request for SLURM jobs
        time_estimate (numeric): time estimate in hours for SLURM jobs
        memory (int): MiXCR memory in GB
    Returns:
        pd.DataFrame: submitted or completed job records

    """
    _validate_backend(backend)

    program_name="MIXCR4.3 Reports"
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
    

    batch_filename = os.path.join(folder, "mixcr_reports_batch.log")
    jobs = []
    for jobname, command in commands.items():
        jobs.append({
            "jobname": jobname,
            "command": command,
            "cwd": folder,
            "log_filename": os.path.join(log_folder, f"{jobname}.log"),
        })

    return _run_mixcr_jobs(
        jobs,
        program_name,
        batch_filename,
        backend=backend,
        max_workers=max_workers,
        cpus=cpus,
        time_estimate=time_estimate,
        memory=memory,
    )


def get_processing_table(folder, show_offtarget=False, offtarget_chain_threshold=0.01):
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
    
    Returns:
        df (pd.DataFrame): dataframe, containing `sample_id`, `extracted_chain` and 
            different processing stats columns. There may be several rows with the same 
            `sample_id`, with each found `extracted_chain`
    """

    if isinstance(folder, list):
        tables = []
        for f in folder:
            table = get_processing_table(f, show_offtarget=show_offtarget)
            tables.append(table)
        return pd.concat(tables).sort_values(by="sample_id").reset_index(drop=True)
    
    results = []
    clonosets = find_all_exported_clonosets_in_folder(folder, chain=None)

    for i, r in clonosets.iterrows():
        sample_id = r["sample_id"]
        chain = r["chain"]
        align_report = read_json_report(sample_id, folder, "align")
        
        try:
            refine_report = read_json_report(sample_id, folder, "refine")
            umi = True
        except FileNotFoundError:
            umi = False
            
        assemble_report = read_json_report(sample_id, folder, "assemble")

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

        results.append([sample_id, chain, Rt, Ru_pc, Ra_pc, Roa_pc, UMIa, UMIc, overseq_threshold, Rf, UMIf, reads_per_umi, Ct, Rcl, Ctc, Rclc, Cfunc, Rfunc, UMIcl, UMIfunc])
    result_df = pd.DataFrame(results, columns=["sample_id", "extracted_chain", "reads_total", "reads_with_umi_pc", "reads_aligned_pc", "reads_overlapped_aln_pc",
                                               "total_umi", "umi_after_correction", "overseq_threshold", "reads_after_filter", "umi_after_filter",
                                               "reads_per_umi", "clones_total", "reads_in_clones_total", "clones", "reads_in_clones", "clones_func", "reads_in_func_clones", "umi_in_clones", "umi_in_func_clones"])
    if not show_offtarget:
        result_df = result_df.loc[result_df.reads_in_clones/result_df.reads_in_clones_total > offtarget_chain_threshold]
    return result_df.sort_values(by="sample_id").reset_index(drop=True)


def show_report_images(folder):
    """
    This function displays QC images `alignQc.svg` and `chainsQc.svg` in Jupyter Notebook.
    This pictures may be generated by `mixcr4_reports` function.
    In case there are no `.svg` images, the `.png` images are shown.

    Args:
        folder (str): folder in which to look for QC images.
    
    Returns:
        None

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


def show_report_images_new(folder, chart_type='summary', count_type='percent', output_file=None):
    """
    Shows quality control reports in MiXCR-like style

    Args:
        folder (str): folder in which to look for QC images.
        chart_type (str): Possible values are `summary` (corresponds to `mixcr exportQc align`, 
            `chains` (`mixcr exportQc chainUsage`)
        count_type (str): possible values are: `percent`, `abs`
        output_file (str): filename ending with '.png' to save an output plot to

    Returns:
        None

    """
    CHAIN_VARIANTS = ['IGH', 'IGK', 'IGL', 'TRA', 'TRB', 'TRD', 'TRL']
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
    for file in files:
        report_type = file[1]
        if report_type == 'align':
            json_report_contents = read_json_report(file[0], folder, report_type=report_type)
            align_data = json_report_contents['notAlignedReasons']
            chain_usage_data = json_report_contents['chainUsage']['chains']
            
            if chart_type == 'summary':
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
                align_df = {}
                for chain, data in chain_usage_data.items(): 
                    align_df.update({chain: data['total'] - data['nonFunctional'],
                                f'{chain} (stops)': data['hasStops'],
                                f'{chain} (OOF)': data['isOOF']})
                df_list.append(pd.DataFrame(align_df, index=[file[0]]))
    results = pd.concat(df_list)
    results = results.sort_index(ascending=False)
    if count_type == 'percent':
        results =  results.div(results.sum(axis=1), axis=0) * 100
    if chart_type == 'summary':
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
        order = sorted(results.columns)
        colors = ['#c26a27', '#ff9429', '#ffcb8f', '#a324b2', '#e553e5', '#faaafa', '#ad3757', '#f05670', '#ffadba', 
                                  '#105bcc', '#2d93fa', '#99ccff', '#198020', '#42b842', '#99e099', '#068a94', '#27c2c2', '#90e0e0', 
                                '#5f31cc', '#845cff', '#c1adff']
        all_variants = ['IGH', 'IGH (stops)', 'IGH (OOF)', 'IGK', 'IGK (stops)', 'IGK (OOF)', 'IGL', 'IGL (stops)', 'IGL (OOF)', 'TRA', 'TRA (stops)', 'TRA (OOF)', 'TRB', 'TRB (stops)', 'TRB (OOF)', 'TRD', 'TRD (stops)', 'TRD (OOF)', 'TRG', 'TRG (stops)', 'TRG (OOF)']
        # colormap = ListedColormap(colors=[colors[i] for i in range(len(colors)) if all_variants[i] in results.columns],
                                    #   name='mixcr')        
        colors = [colors[i] for i in range(len(colors)) if all_variants[i] in results.columns]
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
    if chart_type == 'summary':
        fig.legend(loc='outside upper center',  title='Alignments rate', ncol=3, frameon=False)
    elif chart_type == 'chains':
        fig.legend(loc='outside upper center',  title='Clonal chain usage', ncol=3, frameon=False)
    # plt.tight_layout()
    sns.despine(left=True, bottom=True)
    plt.show()
    if output_file is not None:
        ax.get_figure().savefig(output_file, bbox_inches='tight')
    return
