from .common_functions import print_progress_bar
from .constants import MIXCR
import pandas as pd
import numpy as np
import os
from .slurm import create_slurm_batch_file, run_slurm_command_from_jupyter
from .io import read_json_report, read_clonoset
from .clonosets import find_all_exported_clonosets_in_folder, filter_nonfunctional_clones
from subprocess import Popen, PIPE
from IPython.display import Image, display, SVG


def mixcr4_analyze_batch(sample_df, output_folder, command_template=None,
                         mixcr_path="mixcr", memory=32, time_estimate=1.5, custom_tag_pattern_column=None):
    
    """
    Function for batch runs of MiXCR software using SLURM.
    For each record in the given `sample_df` this function creates a SLURM-script in
    `~/temp/slurm` folder and adds them to SLURM-queue. All the `stdout` logs are also 
    put to `~/temp/slurm` folder. In case of troubles check the latest logs in this folder. 
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

    Returns:
        None
    """
    max_memory = 1500
    min_memory = 16

    program_name="MIXCR4 Analyze Batch"
    samples_num = sample_df.shape[0]

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
    os.makedirs(output_folder, exist_ok=True)

    # default mixcr analyze slurm parameters. They are quite excessive, works fine.
    # time_estimate=1.5
    cpus=40
    if not isinstance(memory, int):
        raise TypeError("memory parameter must be an integer")
    if memory < min_memory:
        print(f"{memory} < than limit ({min_memory}), using {min_memory} GB")
        memory = min_memory
    if memory > max_memory:
        print(f"{memory} > than limit ({max_memory}), using {max_memory} GB")
        memory = max_memory
        
    
    # create slurm batch file for progress tracking
    slurm_batch_filename = os.path.join(output_folder, "mixcr_analyze_slurm_batch.log")
    create_slurm_batch_file(slurm_batch_filename, program_name, samples_num)
    
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
        command = f"cd {output_folder}; " + command
        jobname = f"mixcr_analyze_{sample_id}"
        
        # for batch task finish tracking:
        command += f'; echo "{jobname} finished" >> {slurm_batch_filename}'
        
        # create slurm script and add job to queue, print stdout of sbatch
        stdout, stderr = run_slurm_command_from_jupyter(command, jobname, cpus, time_estimate, memory)
        print(stdout, stderr)
    
    print(f"{samples_num} tasks added to slurm queue\n")
    print(f'To see running progress bar run this function in the next jupyter cell:\nslurm.check_slurm_progress("{slurm_batch_filename}", loop=True)')
    print(f'To see current progress:\nslurm.check_slurm_progress("{slurm_batch_filename}")')


def mixcr_7genes_run_batch(sample_df, output_folder, mixcr_path="mixcr", memory=32, time_estimate=1.5):
    '''
    
    '''
    # default mixcr analyze slurm parameters. They are quite excessive, works fine.
    max_memory = 1500
    min_memory = 16
    cpus=40
    
    program_name="MIXCR4 Analyze 7genes Batch"
    samples_num = sample_df.shape[0]
        
    # Create output dir if does not exist
    os.makedirs(output_folder, exist_ok=True)

    
    if not isinstance(memory, int):
        raise TypeError("memory parameter must be an integer")
    if memory < min_memory:
        print(f"{memory} < than limit ({min_memory}), using {min_memory} GB")
        memory = min_memory
    if memory > max_memory:
        print(f"{memory} > than limit ({max_memory}), using {max_memory} GB")
        memory = max_memory
        
    # create slurm batch file for progress tracking
    slurm_batch_filename = os.path.join(output_folder, "mixcr_analyze_slurm_batch.log")
    create_slurm_batch_file(slurm_batch_filename, program_name, samples_num)
    
    list_of_incomplete_rearrangements = ["DJ_TRB", "VDD_TRD", "DDJ_TRD", "DD_TRD", "DJ_IGH", "VKDE_IGK", "CINTRON_KDE_IGK"]

    # main cycle by samples
    for i,r in sample_df.iterrows():
        sample_id = r["sample_id"]
        r1 = r["R1"]
        r2 = r["R2"]
        output_prefix = sample_id
        
        R1na = f"{sample_id}_R1_not_aligned.fastq.gz"
        R2na = f"{sample_id}_R2_not_aligned.fastq.gz"
        
        commands = [f"cd {output_folder}"]
        
        first_command = f'{mixcr_path} -Xmx{memory}g analyze milab-human-dna-xcr-7genes-multiplex -f --not-aligned-R1 {R1na} --not-aligned-R2 {R2na} {r1} {r2} {output_prefix}'
        commands.append(first_command)
        
        for rearrangement in list_of_incomplete_rearrangements:
            
            # swap r and Rna so we would not implement copy of R_na
            r1, R1na = R1na, r1
            r2, R2na = R2na, r2
            
            output_prefix = os.path.join(rearrangement, sample_id)
            
            R1na = f"{output_prefix}_R1_not_aligned.fastq.gz"
            R2na = f"{output_prefix}_R2_not_aligned.fastq.gz"
            
            i_r_command = f'{mixcr_path} -Xmx{memory}g analyze generic-amplicon -f --species hs --library {rearrangement} --dna --floating-left-alignment-boundary --floating-right-alignment-boundary J -MexportClones.splitFilesBy=[] --not-aligned-R1 {R1na} --not-aligned-R2 {R2na} {r1} {r2} {output_prefix}'
            commands.append(i_r_command)
        
        jobname = f"mixcr_analyze_{sample_id}"
        
        # for batch task finish tracking:
        commands.append(f'echo "{jobname} finished" >> {slurm_batch_filename}')
        
        # join commands by && so that next command runs if previous was finished without error and add new lines to the script
        command = " && \\ \n".join(commands)
        
        # create slurm script and add job to queue, print stdout of sbatch
        stdout, stderr = run_slurm_command_from_jupyter(command, jobname, cpus, time_estimate, memory)
        print(stdout, stderr)
        # print(command)
    print(f"{samples_num} tasks added to slurm queue\n")
    print(f'To see running progress bar run this function in the next jupyter cell:\nslurm.check_slurm_progress("{slurm_batch_filename}", loop=True)')
    print(f'To see current progress:\nslurm.check_slurm_progress("{slurm_batch_filename}")')


def mixcr4_reports(folder, mixcr_path="mixcr"):
    
    """
    runs `mixcr exportQc` commands - `align`, `chainUsage` and `tags` in a given folder 
    for all `.clns` filenames. `align` and `chainUsage` are run twice to create both 
    `svg` and `pdf` files.

    Args:
        folder (str): folder in which to run the `mixcr exportQc` commands
        mixcr_path (str): path to MiXCR binary
    Returns:
        None

    """


    program_name="MIXCR4.3 Reports"
    time_estimate=1
    cpus=40
    memory=32
    
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
    

    
    commands = {"alignQc": f"cd {folder}; {mixcr_path} -Xmx32g exportQc align -f {clns_filenames} {align_filename}",
                "chainUsage": f"cd {folder}; {mixcr_path} -Xmx32g exportQc chainUsage -f {clns_filenames} {chains_filename}",
                "alignQcPDF": f"cd {folder}; {mixcr_path} -Xmx32g exportQc align -f {clns_filenames} {align_filename_pdf}",
                "chainUsagePDF": f"cd {folder}; {mixcr_path} -Xmx32g exportQc chainUsage -f {clns_filenames} {chains_filename_pdf}",
                "tagsQc": f"cd {folder}; {mixcr_path} -Xmx32g exportQc tags -f {clns_filenames} {tags_filename}"#,
                #"postanalysis": f"{MIXCR} -Xmx32g postanalysis individual -f --default-downsampling none --default-weight-function umi --only-productive --tables {tables_filename} --preproc-tables {preproc_filename} {clns_filenames} {postanalysis_filename}"
               }
    

    commands_num = len(commands)
    
    slurm_batch_filename = os.path.join(folder, "mixcr_reports_slurm_batch.log")
    create_slurm_batch_file(slurm_batch_filename, program_name, commands_num)
    
    for jobname, command in commands.items():
        command += f'; echo "{jobname} finished" >> {slurm_batch_filename}'
        stdout, stderr = run_slurm_command_from_jupyter(command, jobname, cpus, time_estimate, memory)
        print(stdout, stderr)


def get_processing_table(folder, show_offtarget=False, off_target_chain_threshold=0.01):
    """
    Searches for clonosets in the the folder, extracts their sample_id's and shows main
    processing stats in a table format. By default does not show "off-target" clonosets - 
    those having less than 1% (default, may be overriden) of reads for the sample_id.
    For example, you have sequenced TRB sample, but there is found 0.5% (by read count) 
    of TRA chains for the same sample_id, then the clonoset will not be shown in the table.
    You can specify `show_offtarget=True` to display all found chains in the table or 
    outherwise set a higher value for `off_target_chain_threshold` (`0.01` by default).

    Args:
        folder (str or list): folder or list of folders in which to look for clonosets and
            processing stats
        show_offtarget (bool): add offtarget chains to the stats
        off_target_chain_threshold (float): threshold for off-target chains
    
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
        clonoset_f = filter_nonfunctional_clones(clonoset)

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
        result_df = result_df.loc[result_df.reads_in_clones/result_df.reads_in_clones_total > off_target_chain_threshold]
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
