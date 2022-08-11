import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
from subprocess import Popen, PIPE
from .common_functions import print_progress_bar
from .constants import MIGEC

def create_barcode_file_from_sample_metadata_subset(input_metadata_df, barcodes_filename):
    metadata_df = input_metadata_df.copy(deep=True)
    metadata_df["Slave barcode sequence"] = ""
    barcode_df = metadata_df[["sample_id", "barcode", "Slave barcode sequence", "raw_data_R1", "raw_data_R2"]].copy(deep=True)
    barcode_df.rename(columns = {"sample_id": "Sample ID", "barcode": "Master barcode sequence",
                                 "raw_data_R1": "", "raw_data_R2": ""}, inplace=True)
    barcode_df.to_csv(barcodes_filename, sep="\t", index=False)
    barcode_df.reset_index(inplace=True, drop=True)
    return barcode_df

def migec_analize(barcodes_filename, output_folder, overseq_threshold=0, suffix="", checkout=True, assemble=True, histogram=True):
    
    # create folder names
    checkout_folder = os.path.join(output_folder, "checkout")
    histogram_folder = os.path.join(output_folder, "histogram")
    assemble_folder = os.path.join(output_folder, "assemble")
    if suffix:
        suffix = "_" + suffix
        checkout_folder += suffix
        histogram_folder += suffix
        assemble_folder += suffix
#     if overseq_threshold:
#         assemble_folder += "_thr_{}".format(overseq_threshold)
    samples_num = len(pd.read_csv(barcodes_filename, sep="\t"))
    
    checkout_log = os.path.join(checkout_folder, "checkout.log")
    histogram_log = os.path.join(histogram_folder, "histogram.log")
    assemble_log = os.path.join(assemble_folder, "assemble.log")
    
    # create output_folder if it doesn't exist
    if not os.path.isdir(output_folder):
        os.mkdir(output_folder)
    
    # run calculations
    print("Running MIGEC Analyze function")
    print("Barcodes filename: {}".format(barcodes_filename))
    print("Samples total: {}".format(samples_num))
    
    
    if checkout:
        migec_checkout(barcodes_filename, output_folder, checkout_folder, checkout_log)
    else:
        print("Checkout process skipped")

    if histogram:
        migec_histogram(samples_num, output_folder, checkout_folder, histogram_folder, histogram_log)
    else:
        print("Histogram process skipped")
    
    if assemble:
        migec_assemble(samples_num, output_folder, checkout_folder, histogram_folder, assemble_folder, overseq_threshold, assemble_log)
    else:
        print("Assemble process skipped")
    
    
def migec_checkout(barcodes_filename, output_folder, checkout_folder, checkout_log=""):
    print("------------------------------")
    print("Running MIGEC Checkout")
    print("Checkout output folder: {}".format(checkout_folder))    
    if not checkout_log:
        checkout_log = os.path.join(checkout_folder, "checkout.log")
    if not os.path.isdir(checkout_folder):
        os.mkdir(checkout_folder)
    samples_num = len(pd.read_csv(barcodes_filename, sep="\t"))
    samples_finished = 0
    
    main_command = "{} CheckoutBatch -cute {} {}".format(MIGEC, barcodes_filename, checkout_folder)
    command = "cd {};".format(output_folder) + main_command
    print("Command: {}".format(main_command))
    
    print_progress_bar(samples_finished, samples_num, program_name="Checkout")
    
    with open(checkout_log, "w") as log_file:
        process = Popen(command, stdout=PIPE, stderr=PIPE, shell=True)
        while True:
            line = process.stdout.readline()
            if not line:
                break
            message = line.decode('ascii')
            if "Finished" in message:
                samples_finished += 1
                print_progress_bar(samples_finished, samples_num, program_name="Checkout")
            log_file.write(message)
        stdout, stderr = process.communicate()
        log_file.write(stderr.decode('ascii'))
        
    print("Saved std output into: {}".format(checkout_log))
    print("MIGEC Checkout Finished")    

def migec_histogram(samples_num, output_folder, checkout_folder, histogram_folder, histogram_log=""):
    print("------------------------------")
    print("Running MIGEC Histogram")
    print("Histogram output folder: {}".format(histogram_folder))
    if not histogram_log:
        histogram_log = os.path.join(histogram_folder, "histogram.log")
    if not os.path.isdir(histogram_folder):
        os.mkdir(histogram_folder)
    main_command = "{} Histogram {} {}".format(MIGEC, checkout_folder, histogram_folder)
    command = "cd {};".format(output_folder) + main_command
    print("Command: {}".format(main_command))
    
    with open(histogram_log, "w") as log_file:
        process = Popen(command, stdout=PIPE, stderr=PIPE, shell=True)
        samples_finished = 0
        
        print_progress_bar(samples_finished, samples_num, program_name="Histogram")
        
        while True:
            line = process.stdout.readline()
            if not line:
                break
            message = line.decode('ascii')
            if "UMIs total" in message:
                samples_finished += 1
                print_progress_bar(samples_finished, samples_num, program_name="Histogram")
            log_file.write(message)
        stdout, stderr = process.communicate()
        log_file.write(stderr.decode('ascii'))
        
    print("Saved std output into: {}".format(histogram_log))
    print("MIGEC Histogram Finished")    
   
    
    
def migec_assemble(samples_num, output_folder, checkout_folder, histogram_folder, assemble_folder, overseq_threshold=0, assemble_log=""):
    print("------------------------------")
    print("Running MIGEC Assemble")
    print("Assemble output folder: {}".format(assemble_folder))
    if not assemble_log:
        assemble_log = os.path.join(assemble_folder, "assemble.log")
    if not os.path.isdir(assemble_folder):
        os.mkdir(assemble_folder)
    overseq_command = ""
    if overseq_threshold:
        print("Using overseq threshold: {}".format(overseq_threshold))
        overseq_command = " --force-overseq {}".format(overseq_threshold)
        
    main_command = "{} AssembleBatch{} {} {} {}".format(MIGEC, overseq_command, checkout_folder, histogram_folder, assemble_folder)
    command = "cd {};" + main_command
    print("Command: {}".format(main_command))
    
    with open(assemble_log, "w") as log_file:
        process = Popen(command, stdout=PIPE, stderr=PIPE, shell=True)
        samples_finished = 0
        
        print_progress_bar(samples_finished, samples_num, program_name="Assemble")
        
        while True:
            line = process.stdout.readline()
            if not line:
                break
            message = line.decode('ascii')
            if "Finished" in message:
                samples_finished += 1
                print_progress_bar(samples_finished, samples_num, program_name="Assemble")
            log_file.write(message)
        stdout, stderr = process.communicate()
        log_file.write(stderr.decode('ascii'))
    
    print("Saved std output into: {}".format(assemble_log))
    print("MIGEC Assemble Finished")

def checkout_stats(folder):
    
    log_filename=os.path.join(folder, "checkout.log.txt")
    checkout_df = pd.read_csv(log_filename, sep="\t")
    files_list = checkout_df["INPUT_FILE_2"].unique()
    sample_ids = []
    good_reads = []
    undef_ms = []
    undef_ss = []
    total_reads_list = []
    percent_good = []
    for file in files_list:
        small_df = checkout_df.loc[checkout_df["INPUT_FILE_2"] == file]
        total_reads = 0
        for i, r in small_df.iterrows():
            value = r["MASTER"]
            name = r["SAMPLE"]
            total_reads += value
            if name == "undef-m":
                undef_ms.append(value)
            elif name == "undef-s":
                undef_ss.append(value)
            else:
                sample_ids.append(name)
                good_reads.append(value)
        total_reads_list.append(total_reads)
        percent_good.append(round(good_reads[-1]/total_reads*100, 2))
    stats_df = pd.DataFrame({"sample_id": sample_ids, "reads_total": total_reads_list, "reads_with_UMI": good_reads,
                             "undef_m": undef_ms, "undef_s": undef_ss, "reads_with_umi_percent": percent_good})
#     stats_df["percent_good"] = (stats_df["reads_with_UMI"])/(stats_df["reads_with_UMI"]+stats_df["undef_m"]+stats_df["undef_s"])*100
    return stats_df
    
def plot_overseq_histograms(histogram_folder):
    df = pd.read_csv(os.path.join(histogram_folder, "overseq.txt"), sep='\t')
    plot_all_overseq_histograms(df, one_plot=True)
    plot_all_overseq_histograms(df, one_plot=False)

    
def plot_all_overseq_histograms(df, one_plot=True):
    if one_plot:
        plot_size = (20,10)
    else:
        plot_size = (10,5)
    plt.rcParams["figure.figsize"] = plot_size
    for indices, row in df.iterrows(): 
        plt.plot(row[2:])
        if not one_plot:
            plt.xlabel(row[0], fontsize=14)
            plt.show()
    if one_plot:
        plt.title("Overseq Histograms", fontsize=18)
        plt.xlabel("Overseq per UMI", fontsize=18)
        plt.show()

