import os
from subprocess import Popen, PIPE
from time import sleep, gmtime, strftime
import re
from .common_functions import print_progress_bar


ALDAN_TEMP_SLURM_DIR = os.path.join(os.path.expanduser("~"), "temp", "slurm")


def partition_by_time(t):
    partition = "short"
    if t >= 2:
        partition = "medium"
    if t >= 16:
        partition = "long"
    if t >= 144:
        partition = "infinite"
    return partition


def run_slurm_command_from_jupyter(command, jobname, cpus, time_estimate, memory, log_filename=None):
    # prepare task parameters and script name
    h=int(time_estimate)
    m=int((time_estimate-h)*60)
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
    partition = partition_by_time(time_estimate)
    
    # create and save slurm script
    slurm_script_template = '''#!/bin/sh

#SBATCH --job-name={jobname}        # Job name
#SBATCH --cpus-per-task={cpus}         # Run on a single CPU
#SBATCH --mem={memory}gb                 # Job memory request
#SBATCH --time={h}:{m}:00           # Time limit hrs:min:sec
#SBATCH --output={log_filename}   # Standard output and error log
#SBATCH --partition={partition}

{command}
    '''
    slurm_text = slurm_script_template.format(jobname=jobname,
                                              cpus=cpus,
                                              memory=memory,
                                              h=h,
                                              m=m,
                                              log_filename=slurm_script_log_filename,
                                              partition=partition,
                                              command=command)  
    
    with open(slurm_script_filename, "w") as f:
        f.write(slurm_text)
        
    # run slurm script
    run_slurm_command = f'chmod a+x {slurm_script_filename}; sbatch {slurm_script_filename}'

    process = Popen(run_slurm_command, stdout=PIPE, stderr=PIPE, shell=True)
    stdout, stderr = process.communicate()
    return stdout, stderr


def create_slurm_batch_file(filename, program_name, tasks_num):
    with open(filename, "w") as f:
        f.write(f"# {program_name}: slurm run {tasks_num} tasks\n")


def check_slurm_progress(filename, loop=False):
    with open(filename, "r") as f:
        first_line = f.readlines()[0]
    match = re.match('# ([^:]+): slurm run ([0-9]+) tasks',first_line)
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