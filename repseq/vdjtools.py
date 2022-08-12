from .common_functions import print_progress_bar
from .constants import VDJTOOLS
import pandas as pd
import os
from subprocess import Popen, PIPE

def vdjtools_convert(folder):
    log_filename = os.path.join(folder, "vdjtools_convert_log.txt")
    metadata_filename = os.path.join(folder, "metadata.txt")
    output_prefix = os.path.join(folder, "vdjtools")
    output_metadata_filename = output_prefix + "_metadata.txt"
    
    metadata_df = pd.read_csv(metadata_filename, sep="\t")
    metadata_df["#file.name"] = output_prefix + "." + metadata_df["sample.id"] + ".txt"
    samples_num = len(metadata_df)
    
    command = "{} Convert -m {} -S MiXcr {}".format(VDJTOOLS, metadata_filename, output_prefix)
      
    print("Running VDJTOOLS Convert")
    print("Convert folder: {}".format(folder))
    print("Command: {}".format(command))
    print("________________________________________________________")
    
    
    with open(log_filename, "w") as log_file:
        process = Popen(command, stdout=PIPE, stderr=PIPE, shell=True)
        samples_finished = 0
        
        print_progress_bar(samples_finished, samples_num, program_name="VDJTOOLS convert")
        
        while True:
            line = process.stdout.readline()
            if not line:
                break
            message = line.decode('ascii')
            if "sample(s).. Writing output" in message:
                samples_finished += 1
                print_progress_bar(samples_finished, samples_num, program_name="VDJTOOLS convert")
            log_file.write(message)
        stdout, stderr = process.communicate()
        log_file.write(stderr.decode('ascii'))
    metadata_df.to_csv(output_metadata_filename, sep="\t", index=False)
    print("________________________________________________________")
    print("VDJTOOLS Convert Finished")    
    print("Saved all metadata for vdjtools to: ", output_metadata_filename)
    print("Written standard output to: {}".format(log_filename))

def filter_non_fuctional_batch(input_folder, output_folder):
    vdjtools_metadata_filename = os.path.join(input_folder, "vdjtools_metadata.txt")
    vdjtools_metadata_raw = pd.read_csv(vdjtools_metadata_filename, sep="\t")
    
    if output_folder == input_folder:
        raise ValueError("input and output folders must be different")
    if not os.path.isdir(output_folder):
        os.mkdir(output_folder)
    
    file_list = []
    for i, r in vdjtools_metadata_raw.iterrows():
        full_filename = r["#file.name"]
        sample_id = r["sample.id"]
        filename = os.path.basename(full_filename)
        new_filename = os.path.join(output_folder, filename)
        clonoset = pd.read_csv(full_filename, sep="\t")
        clonoset = clonoset.loc[~clonoset["cdr3aa"].str.contains("\*|_")]
        clonoset.sort_values(by="count", ascending=False, inplace=True) 
        clonoset.loc[:,"frac"]=clonoset["count"]/clonoset["count"].sum()
        clonoset.to_csv(new_filename, sep="\t", index=False)
        file_list.append((new_filename, sample_id))
    new_vdjtools_metadata = pd.DataFrame(file_list, columns=["#file.name", "sample.id"])
    new_vdjtools_metadata_filename = os.path.join(output_folder, "vdjtools_metadata.txt")
    new_vdjtools_metadata.to_csv(new_vdjtools_metadata_filename, sep="\t", index=False)
    print(f"Written only functional clonosets to {output_folder}")
    print(f"New VDJtools metadata filename: {new_vdjtools_metadata_filename}")