from .common_functions import print_progress_bar
from .constants import MIXCR
import pandas as pd
import numpy as np
import os
from .slurm import create_slurm_batch_file, run_slurm_command_from_jupyter
from .io import read_json_report, read_mixcr_clonoset
from .clonosets import find_all_exported_clonosets_in_folder, filter_nonfunctional_clones
from subprocess import Popen, PIPE
from IPython.display import Image, display, SVG



def mixcr4_analyze_batch(sample_df, output_folder, command_template=None,
                         mixcr_path="mixcr", memory=32, time_estimate=1.5):
    
    """
    Function for batch runs of MiXCR software using SLURM.

    Args:
        sample_df (pd.DataFrame): DataFrame, containing 'sample_id' column and 
            'R1' and 'R2' columns, containing paths (recommended full paths) to raw read files
        output_folder (str): path to output folder
        command_template (str): MiXCR command template 
            (default: 'mixcr analyze milab-human-tcr-rna-multiplex-cdr3 -f r1 r2 output_prefix').
            May be used as an example
        mixcr_path (str): path to MiXCR binary
        memory (int): required OOM in GB
        time_estimate (numeric): time estimate in hours for the calculation. It
            is the limit for SLURM task
    """
    max_memory = 1500
    min_memory = 16

    program_name="MIXCR4.3 Analyze Batch"
    samples_num = sample_df.shape[0]
    
    # by default use the most popular preset for MiLaboratory Human TCR UMI MULTIPLEX Kit
    default_command_template = "mixcr analyze milab-human-tcr-rna-multiplex-cdr3 -f r1 r2 output_prefix"
    if command_template is None:
        command_template = default_command_template
        
    # cut placeholders from command template
    remove_list = ["mixcr", "r1", "r2", "output_prefix"]
    command_template = ' '.join([w for w in command_template.split() if w not in remove_list])
    
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


def mixcr4_reports(folder, mixcr_path="mixcr"):
    
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
        clonoset = read_mixcr_clonoset(r.filename)
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
            overseq_threshold = int(refine_report["correctionReport"]["filterReport"]["operatorReports"][0]["operatorReport"]["threshold"])
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

        results.append([sample_id, chain, Rt, Ru_pc, Ra_pc, Roa_pc, UMIa, UMIc, overseq_threshold, Rf, UMIf, reads_per_umi, Ct, Rcl, Ctc, Rclc, Cfunc, Rfunc, UMIcl, UMIfunc])
    result_df = pd.DataFrame(results, columns=["sample_id", "extracted_chain", "reads_total", "reads_with_umi_pc", "reads_aligned_pc", "reads_overlapped_aln_pc",
                                               "total_umi", "umi_after_correction", "overseq_threshold", "reads_after_filter", "umi_after_filter",
                                               "reads_per_umi", "clones_total", "reads_in_clones_total", "clones", "reads_in_clones", "clones_func", "reads_in_func_clones", "umi_in_clones", "umi_in_func_clones"])
    if not show_offtarget:
        result_df = result_df.loc[result_df.reads_in_clones/result_df.reads_in_clones_total > off_target_chain_threshold]
    return result_df.sort_values(by="sample_id").reset_index(drop=True)

def show_report_images(folder):
    """
    This function displays 

    Args:
        folder (str): folder in which to look for QC images. 
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
        