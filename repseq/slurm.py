import math
import os
import re
import shutil
from subprocess import Popen, PIPE
from time import sleep, gmtime, strftime
from .common_functions import print_progress_bar


ALDAN_TEMP_SLURM_DIR = os.path.join(os.path.expanduser("~"), "temp", "slurm")
STANDARD_PARTITIONS = ("short", "medium", "long", "infinite")
PARTITION_LIMIT_HOURS = {}
SLURM_NOT_AVAILABLE_MESSAGE = "SLURM is not installed or 'sinfo' is not available on PATH."
SINFO_FORMAT = "%P|%l|%N|%f"


def _run_sinfo():
    if shutil.which("sinfo") is None:
        raise RuntimeError(SLURM_NOT_AVAILABLE_MESSAGE)

    process = Popen(
        ["sinfo", "--all", "--Node", "--noheader", f"--format={SINFO_FORMAT}"],
        stdout=PIPE,
        stderr=PIPE,
    )
    stdout, stderr = process.communicate()
    if process.returncode != 0:
        error_text = stderr.decode(errors="replace").strip()
        message = "Unable to query SLURM partitions with sinfo."
        if error_text:
            message += f" {error_text}"
        raise RuntimeError(message)
    return stdout.decode(errors="replace")


def _parse_sinfo_output(output):
    rows = []
    for line in output.splitlines():
        if not line.strip():
            continue
        parts = line.split("|", 3)
        if len(parts) != 4:
            continue
        partition, time_limit, nodes, constraints = (part.strip() for part in parts)
        rows.append({
            "partition": partition.rstrip("*"),
            "default": partition.endswith("*"),
            "time_limit": time_limit,
            "nodes": nodes,
            "constraints": constraints,
        })
    return rows


def _slurm_time_to_hours(time_limit):
    normalized = str(time_limit).strip().lower()
    if normalized in {"infinite", "unlimited"}:
        return math.inf
    if not normalized:
        raise ValueError("Empty SLURM time limit")

    days = 0
    time_text = normalized
    if "-" in normalized:
        days_text, time_text = normalized.split("-", 1)
        days = int(days_text)

    fields = time_text.split(":")
    if len(fields) == 3:
        hours, minutes, seconds = fields
    elif len(fields) == 2:
        hours = 0
        minutes, seconds = fields
    elif len(fields) == 1:
        hours = 0
        minutes = fields[0]
        seconds = 0
    else:
        raise ValueError(f"Unsupported SLURM time limit: {time_limit}")

    return days * 24 + int(hours) + int(minutes) / 60 + int(seconds) / 3600


def _get_sinfo_rows():
    return _parse_sinfo_output(_run_sinfo())


def update_partition_limits():
    """Query sinfo and refresh limits for short, medium, long, and infinite."""
    discovered_limits = {}
    for row in _get_sinfo_rows():
        partition = row["partition"]
        if partition not in STANDARD_PARTITIONS:
            continue
        try:
            limit_hours = _slurm_time_to_hours(row["time_limit"])
        except (TypeError, ValueError):
            continue
        previous_limit = discovered_limits.get(partition)
        if previous_limit is None or limit_hours > previous_limit:
            discovered_limits[partition] = limit_hours

    if not discovered_limits:
        names = ", ".join(STANDARD_PARTITIONS)
        raise RuntimeError(f"None of the standard SLURM partitions ({names}) were found by sinfo.")

    PARTITION_LIMIT_HOURS.clear()
    PARTITION_LIMIT_HOURS.update(discovered_limits)
    return PARTITION_LIMIT_HOURS.copy()


def partition_by_time(time_estimate, partition_limits=None):
    """Return the smallest discovered standard partition that fits the requested hours."""
    try:
        requested_hours = float(time_estimate)
    except (TypeError, ValueError):
        raise TypeError("time_estimate must be a positive number of hours") from None
    if not math.isfinite(requested_hours) or requested_hours <= 0:
        raise ValueError("time_estimate must be a positive finite number of hours")

    if partition_limits is None:
        partition_limits = update_partition_limits()

    candidates = []
    for order, partition in enumerate(STANDARD_PARTITIONS):
        if partition not in partition_limits:
            continue
        limit = partition_limits[partition]
        if isinstance(limit, str):
            limit = _slurm_time_to_hours(limit)
        else:
            limit = float(limit)
        if requested_hours <= limit:
            candidates.append((limit, order, partition))

    if not candidates:
        limits_text = ", ".join(
            f"{partition}={partition_limits[partition]}"
            for partition in STANDARD_PARTITIONS
            if partition in partition_limits
        )
        raise ValueError(
            f"Requested time ({requested_hours:g} hours) exceeds available SLURM partition limits"
            + (f": {limits_text}" if limits_text else ".")
        )
    return min(candidates)[2]


def _constraint_names(rows):
    names = set()
    ignored_values = {"", "(null)", "none", "n/a"}
    for row in rows:
        for name in row["constraints"].split(","):
            name = name.strip()
            if name.lower() not in ignored_values:
                names.add(name)
    return sorted(names)


def info():
    """Print SLURM partitions, limits, nodes, and available constraint names."""
    try:
        rows = _get_sinfo_rows()
    except RuntimeError as error:
        print(error)
        return

    if not rows:
        print("SLURM sinfo returned no partition or node information.")
        return

    grouped = {}
    for row in rows:
        partition = row["partition"]
        details = grouped.setdefault(partition, {
            "default": False,
            "limits": set(),
            "nodes": set(),
            "constraints": set(),
        })
        details["default"] = details["default"] or row["default"]
        if row["time_limit"]:
            details["limits"].add(row["time_limit"])
        if row["nodes"]:
            details["nodes"].add(row["nodes"])
        if row["constraints"].lower() not in {"", "(null)", "none", "n/a"}:
            details["constraints"].add(row["constraints"])

    headers = ("Partition", "Default", "Time limit", "Nodes", "Constraints")
    table_rows = []
    for partition, details in grouped.items():
        table_rows.append((
            partition,
            "yes" if details["default"] else "",
            ", ".join(sorted(details["limits"])),
            ", ".join(sorted(details["nodes"])),
            ", ".join(sorted(details["constraints"])),
        ))

    widths = [
        max(len(headers[index]), *(len(row[index]) for row in table_rows))
        for index in range(len(headers))
    ]
    print("SLURM partitions and nodes:")
    print("  ".join(header.ljust(widths[index]) for index, header in enumerate(headers)))
    print("  ".join("-" * width for width in widths))
    for row in table_rows:
        print("  ".join(value.ljust(widths[index]) for index, value in enumerate(row)))

    constraints = _constraint_names(rows)
    if constraints:
        print(f"Possible constraints: {', '.join(constraints)}")
    else:
        print("Possible constraints: none reported by sinfo")


def run_slurm_command_from_jupyter(command, jobname=None, cpus=40, time_estimate=1.5,
                                    memory=32, log_filename=None, verbose=True,
                                    constraint=None):
    """Create and submit a SLURM script using live sinfo partition limits.

    Args:
        command (str): Shell command to run in the submitted job.
        jobname (str): Short SLURM job name.
        cpus (int): CPUs requested for the task.
        time_estimate (numeric): Requested runtime in hours.
        memory (int): Requested memory in GB.
        log_filename (str): Optional output log filename in an existing directory.
        verbose (bool): Print the submission output when true.
        constraint (str): Optional SLURM constraint expression, such as ``hpc``.

    Returns:
        tuple: ``stdout`` and ``stderr`` returned by ``sbatch``.
    """
    # prepare task parameters and script name
    if not isinstance(jobname, str) or jobname == "":
        raise ValueError("Please set a short meaningful jobname")
    if constraint is not None:
        if not isinstance(constraint, str) or not constraint.strip():
            raise ValueError("constraint must be a non-empty string or None")
        if "\n" in constraint or "\r" in constraint:
            raise ValueError("constraint must not contain line breaks")

    partition = partition_by_time(time_estimate)
    requested_hours = float(time_estimate)
    total_seconds = math.ceil(requested_hours * 3600)
    h, remainder = divmod(total_seconds, 3600)
    m, s = divmod(remainder, 60)

    os.makedirs(ALDAN_TEMP_SLURM_DIR, exist_ok=True)
    datetime = strftime("%Y_%m_%d__%H_%M_%S", gmtime())
    slurm_script_filename = os.path.join(ALDAN_TEMP_SLURM_DIR, f"{datetime}_{jobname}.sh")
    if isinstance(log_filename, str):
        log_dirname = os.path.dirname(log_filename)
        if os.path.isdir(log_dirname):
            slurm_script_log_filename = log_filename
        else:
            raise FileNotFoundError(f"dir '{log_dirname}' not found for log_file")
    else:
        slurm_script_log_filename = os.path.join(ALDAN_TEMP_SLURM_DIR, f"{datetime}_{jobname}.log")

    constraint_line = "" if constraint is None else f"#SBATCH --constraint={constraint}\n"

    # create and save slurm script
    slurm_script_template = '''#!/bin/sh

#SBATCH --job-name={jobname}        # Job name
#SBATCH --cpus-per-task={cpus}         # Run on a single CPU
#SBATCH --mem={memory}gb                 # Job memory request
#SBATCH --time={h}:{m:02d}:{s:02d}           # Time limit hrs:min:sec
#SBATCH --output={log_filename}   # Standard output and error log
#SBATCH --partition={partition}
{constraint_line}
{command}
    '''
    slurm_text = slurm_script_template.format(jobname=jobname,
                                              cpus=cpus,
                                              memory=memory,
                                              h=h,
                                              m=m,
                                              s=s,
                                              log_filename=slurm_script_log_filename,
                                              partition=partition,
                                              constraint_line=constraint_line,
                                              command=command)

    with open(slurm_script_filename, "w") as f:
        f.write(slurm_text)

    # run slurm script
    run_slurm_command = f'chmod a+x {slurm_script_filename} && sbatch {slurm_script_filename}'
    process = Popen(run_slurm_command, stdout=PIPE, stderr=PIPE, shell=True)
    stdout, stderr = process.communicate()
    if verbose:
        print(jobname, stdout, stderr)
    return stdout, stderr


def create_slurm_batch_file(filename, program_name, tasks_num):
    with open(filename, "w") as f:
        f.write(f"# {program_name}: slurm run {tasks_num} tasks\n")


def check_slurm_progress(filename, loop=False):
    with open(filename, "r") as f:
        first_line = f.readlines()[0]
    match = re.match('# ([^:]+): slurm run ([0-9]+) tasks', first_line)
    program_name = match[1]
    tasks = int(match[2])

    if loop:
        while True:
            with open(filename, "r") as f:
                text = f.read()
            finished = len(re.findall("finished", text))
            sleep(0.5)
            print_progress_bar(finished, tasks, program_name=program_name, object_name="task(s)")
            if finished == tasks:
                break

    else:
        with open(filename, "r") as f:
            text = f.read()
        finished = len(re.findall("finished", text))
        print(text)
        print_progress_bar(finished, tasks, program_name=program_name, object_name="task(s)")
