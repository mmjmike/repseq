from .common_functions import print_progress_bar
from .constants import MIXCR
import pandas as pd
import numpy as np
import os
from subprocess import Popen, PIPE

def mixcr_analyze(assemble_results_dir, sample_metadata=[], chain="TRB"):
    #parameters of mixcr analyze are here https://mixcr.readthedocs.io/en/master/analyze.html
    

    with open(metadata_filename,"w") as meta_file, open(log_filename,"w") as log_file:
        
        for index, row in assembled_samples_df.iterrows():
            sample=row["#SAMPLE_ID"]

            if not chain_fixed:
                if sample not in list(sample_metadata["sample_id"]):
                    continue
                chain = sample_metadata.loc[sample_metadata["sample_id"] == sample]["chain"].iloc[0]
                receptor_type = "--receptor-type {}".format(chain.upper())
            #if not sample.startswith("UCB"): #some filter for sample names can be here
            #    continue
            R1=row["OUTPUT_ASSEMBLY1"]
            R2=row["OUTPUT_ASSEMBLY2"]
            log_file.write("Analyzing sample: {}\nCommand:\n".format(sample))
            sample_filename_prefix = os.path.join(assemble_results_dir, sample)
            
            command="""{} analyze amplicon -f \
                            --species hs \
                            --starting-material rna \
                            --5-end no-v-primers --3-end c-primers \
                             {} \
                            --adapters no-adapters \
                            --report {}_report.txt {} {} {}""".format(MIXCR, receptor_type,
                           sample_filename_prefix,
                           R1,
                           R2,
                           sample_filename_prefix)

            log_file.write(command)

            process = Popen(command, stdout=PIPE, stderr=PIPE, shell=True)
            
            stdout, stderr = process.communicate()
            stdout_message = stdout.decode('ascii')
            stderr_message = stderr.decode('ascii')
            log_file.write(stdout_message)
            log_file.write(stderr_message)
            
            meta_file_line = "\n{}.clonotypes.{}.txt\t{}".format(sample, chain.upper(), sample)
            meta_file.write(meta_file_line)
            
            print_progress_bar(samples_finished, samples_num, program_name="MIXCR Analyze Batch")
    print("Written metadata to: {}".format(metadata_filename))
    print("Written standard output to: {}".format(log_filename))

    
