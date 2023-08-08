from .common_functions import print_progress_bar
from .constants import MIXCR
import pandas as pd
import numpy as np
import os
from .slurm import create_slurm_batch_file, run_slurm_command_from_jupyter
from .io import read_json_report, read_mixcr_clonoset
from .clonosets import find_all_exported_clonosets_in_folder, filter_nonfunctional_clones
from subprocess import Popen, PIPE
from IPython.display import Image, display



def mixcr4_analyze_batch(sample_df, output_folder, command_template=None, mixcr_path="mixcr"):
    
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
    time_estimate=1.5
    cpus=40
    memory=32
    
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
    align_filename = "alignQc.png"
    chains_filename = "chainsQc.png"
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


def get_processing_table(folder, show_offtarget=False):
    
    if isinstance(folder, list):
        tables = []
        for f in folder:
            table = get_processing_table(f, show_offtarget=show_offtarget)
            tables.append(table)
        return pd.concat(tables).sort_values(by="sample_id").reset_index(drop=True)
    
    results = []
    clonosets = find_all_exported_clonosets_in_folder(folder, chain=None)
    
    off_target_chain_threshold = 0.01

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
    images = [os.path.join(folder, "alignQc.png"),
          os.path.join(folder, "chainsQc.png")]
    for imageName in images:
        print(imageName)
        display(Image(filename=imageName))


################ Old functions
def mixcr_check_input(input_folder, output_folder, protocol, samples_metadata_df, chain, species):
    os.makedirs(output_folder, exist_ok=True)
    
    protocols = ["RACE", "multiplex", "umi_multiplex", "shotgun"]
    chains = ["tra", "trb", "trg", "trd", "tcr", "bcr", "igh", "igl", "igk"]
    species_variants = ["hs", "mmu", "hsa", "musmusculus", "HomoSapiens"]
    
    # Check protocol
    if protocol not in protocols:
        error_msg = "Unknown protocol '{}'. Possible values: {}".format(protocol, ", ".join(protocols))
        return False, error_msg, None
    assembly = True
    if protocol in ["multiplex", "shotgun"]:
        assembly = False
        
    # Check species
    species = species.lower()
    if species not in species_variants:
        return False, "ERROR. Unknown species '{}'. Choose one of the following: {}".format(species, " ".join(species_variants)), None
         
    # Check presence of sample dataframe if the raw data is not assembled with MiGEC or MiNNN
    if not assembly and samples_metadata_df is None:
        return False, "ERROR! If the data is not assembly with MiGEC or MiNNN, 'samples_metadata_df' parameter is needed", None
    
    # Check chains
    chain_from_df = False
    if isinstance(samples_metadata_df, pd.DataFrame):
        samples_metadata_df.rename(columns={"Chain": "chain",
                                            "R1": "raw_data_R1",
                                            "R2": "raw_data_R2",
                                            "r1": "raw_data_R1",
                                            "r2": "raw_data_R2"}, inplace=True)
        if "chain" in samples_metadata_df.columns:
            chain_from_df = True
            samples_metadata_df[["chain"]] = samples_metadata_df[["chain"]].fillna(value="")
            for df_chain in samples_metadata_df.chain:
                if df_chain.lower() not in chains:
                    return False, "ERROR. Unknown chain in dataframe '{}'. Choose one of the following: {}".format(df_chain, " ".join(chains)), None
        if not assembly:
            if "raw_data_R1" not in samples_metadata_df.columns:
                return False, "ERROR. R1 raw data is not specified in dataframe", None
            if "raw_data_R2" not in samples_metadata_df.columns:
                return False, "ERROR. R1 raw data is not specified in dataframe", None             
    if chain is None:
        if not chain_from_df:
            return False, "ERROR. Chain is not specified neither in 'chain' parameter, nor in dataframe", None
    elif chain.lower() not in chains:
        return False, "ERROR. Unknown chain '{}'. Choose one of the following: {}".format(chain, " ".join(chains)), None
    
    # If assembly - read "assemble.log.txt" from input folder, and specify chain in samples_metadata
    if assembly:
        if not os.path.isdir(input_folder):
            return False, f"ERROR. Input folder does not exist: '{input_folder}'", None
        input_filename = os.path.join(input_folder, "assemble.log.txt")
        if not os.path.isfile(input_filename):
            return False, f"ERROR. Input file does not exist: '{input_folder}'", None
        samples_metadata=pd.read_csv(input_filename,
                               sep="\t",header=0).rename(columns={"#SAMPLE_ID": "sample_id",
                                                                  "OUTPUT_ASSEMBLY1": "raw_data_R1",
                                                                  "OUTPUT_ASSEMBLY2": "raw_data_R2"})
        if chain_from_df:
            samples_metadata = samples_metadata.merge(samples_metadata_df[["sample_id", "chain"]])
        else:
            samples_metadata["chain"] = chain.lower()
    else:
        samples_metadata = samples_metadata_df
    # check R1 and R2 filenames. Add full paths, if there are local paths
    
    # try find at least one global path:
    global_path_exists = False
    for i, r in samples_metadata.iterrows():
        r1 = r["raw_data_R1"]
        r2 = r["raw_data_R2"]
        if os.path.isfile(r1) or os.path.isfile(r2):
            global_path_exists = True
            break
    # try to add global path, if there were none
    if not global_path_exists:
        if not os.path.isdir(input_folder):
            return False, f"ERROR. Input folder does not exist: '{input_folder}'", None
        samples_metadata["raw_data_R1"] = samples_metadata["raw_data_R1"].apply(lambda x: os.path.join(input_folder, x))
        samples_metadata["raw_data_R2"] = samples_metadata["raw_data_R2"].apply(lambda x: os.path.join(input_folder, x))
    
    # check presence of all raw R1/R2 input files
    for i, r in samples_metadata.iterrows():
        r1 = r["raw_data_R1"]
        r2 = r["raw_data_R2"]
        if not os.path.isfile(r1):
            return False, f"ERROR. Raw input file does not exist: '{r1}'", None
        if not os.path.isfile(r2):
            return False, f"ERROR. Raw input file does not exist: '{r2}'", None
    return True, "", samples_metadata

def mixcr_create_common_command(protocol, species, add_params):
    
    protocol_parameters = {"RACE":          {"analyze_mode": "amplicon", "5end": "--5-end no-v-primers", "3end": "--3-end c-primers",  "adapters": "--adapters no-adapters"},
                           "multiplex":     {"analyze_mode": "amplicon", "5end": "--5-end v-primers",  "3end": "--3-end c-primers",  "adapters": "--adapters adapters-present"},
                           "umi_multiplex": {"analyze_mode": "amplicon", "5end": "--5-end v-primers",  "3end": "--3-end c-primers",  "adapters": "--adapters no-adapters"},
                           "shotgun":       {"analyze_mode": "shotgun",  "5end": "", "3end": "", "adapters": ""}}
    
    mixcr_path = "/software/bin/mixcr"
    
    additional_parameters = []
    if add_params["contig_assembly"]:
        additional_parameters.append("--contig-assembly")
    if add_params["only_productive"]:
        additional_parameters.append("--only-productive")
    if add_params["impute_germline"]:
        additional_parameters.append("--impute-germline-on-export")
    if add_params["save_reads"]:
        additional_parameters.append('--"-OsaveOriginalReads=true"')
    
    command = """{mixcr} analyze {amplicon_or_shotgun} -f --species {species} --starting-material rna \
    {primers_5} {primers_3} --receptor-type {{chain}} {adapters} {add_parameters} --report {{output_prefix}}_report.txt \
    {{r1}} {{r2}} {{output_prefix}}""".format(mixcr=MIXCR, 
                                                                      amplicon_or_shotgun=protocol_parameters[protocol]["analyze_mode"],
                                                                      species=species,
                                                                      primers_5=protocol_parameters[protocol]["5end"],
                                                                      primers_3=protocol_parameters[protocol]["3end"],
                                                                      adapters=protocol_parameters[protocol]["adapters"],
                                                                      add_parameters=" ".join(additional_parameters))
    return command
    
def mixcr_run_batch(common_command, samples_metadata, input_folder, output_folder):
    program_name="MIXCR Analyze Batch"
    
    bcr_variants = ["igh", "igl", "igk"]
    tcr_variants = ["tra", "trb", "trg", "trd"]
    
    samples_num = samples_metadata.shape[0]
    start_message = "Starting MiXCR batch analysys.\nTotal samples:\t{}\nCommon command: {}".format(str(samples_num), common_command)
    print(start_message)
    print("_________________________________________________________________")

    samples_finished = 0
    print_progress_bar(samples_finished, samples_num, program_name=program_name)
    output_files = []
    log_filename = os.path.join(output_folder, "mixcr_analyze.log")
    
    with open(log_filename,"w") as log_file:

        for index, row in samples_metadata.iterrows():
            sample=row["sample_id"]
            chain = row["chain"]
            log_file.write("Analyzing sample: {}\n".format(sample))
            output_prefix = os.path.join(output_folder, sample)
            command=common_command.format(chain=chain,
                                          output_prefix=output_prefix,
                                          r1=row["raw_data_R1"],
                                          r2=row["raw_data_R2"])
            
            if chain.lower() == "bcr":
                for bcr_variant in bcr_variants:
                    filename = "{}.clonotypes.{}.txt".format(sample, bcr_variant.upper())
                    new_sample_name = "{}_{}".format(sample, bcr_variant.upper())
                    output_files.append((filename, new_sample_name, sample))
            elif chain.lower() == "tcr":
                for tcr_variant in tcr_variants:
                    filename = "{}.clonotypes.{}.txt".format(sample, tcr_variant.upper())
                    new_sample_name = "{}_{}".format(sample, tcr_variant.upper())
                    output_files.append((filename, new_sample_name, sample))
            else:
                filename = "{}.clonotypes.{}.txt".format(sample, chain.upper())
                output_files.append((filename, sample))
    
            log_file.write(command)
            process = Popen(command, stdout=PIPE, stderr=PIPE, shell=True)
            stdout, stderr = process.communicate()
            stdout_message = stdout.decode('ascii')
            stderr_message = stderr.decode('ascii')
            log_file.write(stdout_message)
            log_file.write(stderr_message)
    
            samples_finished += 1
            print_progress_bar(samples_finished, samples_num, program_name=program_name)
    print("Written standard output to: {}".format(log_filename))
    
    output_files = [sample_data for sample_data in output_files if not clonoset_is_empty(os.path.join(output_folder, sample_data[0]))]
    return output_files

def clonoset_is_empty(filename):
    if not os.path.isfile(filename):
        return False
    df= pd.read_csv(filename, sep="\t")
    return len(df) == 0

def mixcr_create_file_list(output_folder, output_files):
    # check if there is at least one multichain sample
    multichain = False
    for file_data in output_files:
        if len(file_data) == 3:
            multichain = True
            break
    columns = ["#file.name", "sample.id"]
    if multichain:
        columns += ["old.sample.id"]
    output_files_df = pd.DataFrame(output_files, columns=columns)
    if multichain:
        output_files_df["old.sample.id"] = output_files_df["old.sample.id"].fillna(output_files_df["sample.id"])
    metadata_filename = os.path.join(output_folder, "metadata.txt")
    output_files_df.to_csv(metadata_filename, sep="\t", index=False)
    print("Written metadata to: {}".format(metadata_filename))

def mixcr_analyze(input_folder, output_folder, protocol, samples_metadata_df=None, chain=None, species="hsa",
                  contig_assembly=False, impute_germline=False, save_reads=False, only_productive=False):
    
    check_input, error_msg, samples_metadata = mixcr_check_input(input_folder, output_folder, protocol, samples_metadata_df, chain, species)
    if not check_input:
        print(error_msg)
        return
    
    additional_params = {"contig_assembly": contig_assembly,
                         "impute_germline": impute_germline,
                         "save_reads": save_reads,
                         "only_productive": only_productive}
    
    command = mixcr_create_common_command(protocol, species, additional_params)
    
    output_files = mixcr_run_batch(command, samples_metadata, input_folder, output_folder)
    
    mixcr_create_file_list(output_folder, output_files)

    
