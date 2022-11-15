def minnn_batch(barcodes_filename, working_dir, overseq_threshold=1, suffix=None):
    barcodes_df = pd.read_csv(barcodes_filename, sep="\t").set_axis(['sample_id', 'barcode', '-', 'R1', 'R2'], axis=1)
    minnn_folder = os.path.join(working_dir, "minnn")
    if not os.path.isdir(minnn_folder):
        os.mkdir(minnn_folder)
    
    minnn_processing(barcodes_df, minnn_folder, overseq_threshold)
    minnn_stats(barcodes_df, minnn_folder, overseq_threshold)
    
def minnn_stats(df, minnn_folder, overseq_threshold):
    print("_"*100)
    print("Calculating processing stats")
    print("_"*100)
    
    stats_df = calculate_reads(df)
    sort_folder = os.path.join(minnn_folder, "sort")
    sort2_folder = os.path.join(minnn_folder, "sort2")
    filter_folder = os.path.join(minnn_folder, "filter")
    sort_stats_folder = os.path.join(minnn_folder, "sort_stats")
    sort2_stats_folder = os.path.join(minnn_folder, "sort2_stats")
    filter_stats_folder = os.path.join(minnn_folder, "filter_stats")
    
    minnn_common("stat-groups", df, input_folder=sort_folder, output_folder=sort_stats_folder)
    minnn_common("stat-groups", df, input_folder=sort2_folder, output_folder=sort2_stats_folder)
    minnn_extract_stats = read_minnn_stats(df, sort_stats_folder).rename(columns={"reads":"reads_extract", "migs":"migs_extract"})
    minnn_correct_stats = read_minnn_stats(df, sort2_stats_folder).rename(columns={"reads":"reads_correct", "migs":"migs_correct"})
    stats_df = stats_df.merge(minnn_extract_stats)
    stats_df = stats_df.merge(minnn_correct_stats)
    stats_df["overseq_threshold"] = overseq_threshold
    
    if overseq_threshold > 1:
        minnn_common("stat-groups", df, input_folder=filter_folder, output_folder=filter_stats_folder)
        minnn_filter_stats = read_minnn_stats(df, filter_stats_folder).rename(columns={"reads":"reads_filter", "migs":"migs_filter"})
        stats_df = stats_df.merge(minnn_filter_stats)
        
    stats_filename = os.path.join(minnn_folder, "minnn_stats.txt")
    stats_df.to_csv(stats_filename, sep="\t", index=False)
    print("Processing stats written to:", stats_filename)
    
    print("_"*100)
    print("Calculating Overseq Histogram")
    print("_"*100)
    minnn_histogram(df, sort2_stats_folder, minnn_folder, reads=True)
    
    return stats_df
    
    
def read_minnn_stats(df, folder):
    results = []
    for i, r in df.iterrows():
        sample_id = r["sample_id"]
        filename=os.path.join(folder, f"{sample_id}.txt")
        sample_stats = pd.read_csv(filename, sep=" ")
        migs = len(sample_stats)
        reads = sample_stats["count"].sum()
        results.append((sample_id, reads, migs))
    return pd.DataFrame(results, columns=["sample_id", "reads", "migs"])

def minnn_histogram(df, input_folder, output_folder, reads=True):
    sample_counts = {}
    df["filename"] = df["sample_id"].apply(lambda x: os.path.join(input_folder, x + ".txt"))
    
    hist_df = overseq_hist_from_df(df, reads=reads)
    
    output_filename = os.path.join(output_folder, "overseq.txt")
    hist_df.fillna(0).to_csv(output_filename, sep="\t", index=False)
    print("Overseq histogram written to:", output_filename)
    
def show_minnn_stats(folder):
    path = os.path.join(folder, "minnn_stats.txt")
    return pd.read_csv(path, sep="\t")
    
def calculate_reads(df):
    results = []
    for i, r in df.iterrows():
        results.append((r["sample_id"], count_reads_in_fastq_gz(r["R1"])))
    return pd.DataFrame(results, columns=["sample_id", "reads_total"])
        
def count_reads_in_fastq_gz(filename):
    with gzip.open(filename, 'rb') as f:
        for i, l in enumerate(f):
            pass
    return int((i+1)/4)
    
def minnn_processing(df, minnn_folder, overseq_threshold):
    # folder names
    extract_folder = os.path.join(minnn_folder, "extract")
    sort_folder = os.path.join(minnn_folder, "sort")
    correct_folder = os.path.join(minnn_folder, "correct")
    sort2_folder = os.path.join(minnn_folder, "sort2")
    filter_folder = os.path.join(minnn_folder, "filter")
    consensus_folder = os.path.join(minnn_folder, "consensus")
    mif2fastq_folder = os.path.join(minnn_folder, "assemble")
    histogram_folder = os.path.join(minnn_folder, "histogram")
    
    # run all commands
    minnn_common("extract", df, input_folder="", output_folder=extract_folder)
    minnn_common("sort", df, input_folder=extract_folder, output_folder=sort_folder)
    minnn_common("correct", df, input_folder=sort_folder, output_folder=correct_folder)
    minnn_common("sort", df, input_folder=correct_folder, output_folder=sort2_folder)
    if overseq_threshold > 1:
        # filter MIGs, that have not enough reads
        minnn_common("filter-by-count", df, input_folder=sort2_folder, output_folder=filter_folder,threshold=overseq_threshold)
        input_concensus_folder = filter_folder
    else:
        input_concensus_folder = sort2_folder
    minnn_common("consensus", df, input_folder=input_concensus_folder, output_folder=consensus_folder)
    minnn_common("mif2fastq", df, input_folder=consensus_folder, output_folder=mif2fastq_folder)
    
def minnn_common(procedure, df, input_folder, output_folder, threshold=1):
    procedures = {"extract": '{minnn} extract --pattern "{barcode}" --input {raw_R1} {raw_R2} --output {output_folder}/{sample_id}.mif',
                "sort": '{minnn} sort --groups UMI --input {input_folder}/{sample_id}.mif --output {output_folder}/{sample_id}.mif',
                "stat-groups": '{minnn} stat-groups --groups UMI --input {input_folder}/{sample_id}.mif --output {output_folder}/{sample_id}.txt',
                "correct": '{minnn} correct --groups UMI --input {input_folder}/{sample_id}.mif --output {output_folder}/{sample_id}.mif',
                "filter-by-count": '{minnn} filter-by-count --groups UMI --min-count {threshold} --input {input_folder}/{sample_id}.mif --output {output_folder}/{sample_id}.mif',
                "consensus": '{minnn} consensus --groups UMI --input {input_folder}/{sample_id}.mif --output {output_folder}/{sample_id}.mif',
                "mif2fastq": '{minnn} mif2fastq --input {input_folder}/{sample_id}.mif --group R1={output_folder}/{sample_id}_R1.fastq R2={output_folder}/{sample_id}_R2.fastq'}
    
    if procedure not in procedures:
        msg = "Unknown MiNNN procedure '{}'. Possible values: {}".format(procedure, ", ".join(procedures))
        raise ValueError(msg)
    
    if not os.path.isdir(output_folder):
        os.mkdir(output_folder)
    
    samples_total = len(df)
    samples_done = 0
    program_name = f"MiNNN {procedure} Batch"
    some_barcode = "^(barcode:tggtatcaacgcagagt(UMI:NtNNNNtNNNNtNNNN))\*"
    common_command = procedures[procedure]
    command_example = common_command.format(minnn=MINNN,
                                            threshold=threshold,
                                            input_folder=input_folder, output_folder=output_folder,
                                            raw_R1="r1_path", raw_R2="r2_path",
                                            barcode=some_barcode, sample_id="sample_id")
    print("_"*100)
    print("Running", program_name)
    print("Command:", command_example)
    
    print_progress_bar(samples_done, samples_total, program_name=program_name)
    
    log_filename = os.path.join(output_folder, f"minnn_{procedure}.log")
    if procedure == "mif2fastq":
        list_of_files = []
    
    with open(log_filename, "w") as log_file:
        for i, r in df.iterrows():
            sample_id = r["sample_id"]
            r1 = r["R1"]
            r2 = r["R2"]
            barcode = r["barcode"]
            log_file.write("Analyzing sample: {}\nCommand:\n".format(sample_id))
            command = common_command.format(minnn=MINNN,
                                            input_folder=input_folder,
                                            output_folder=output_folder,
                                            raw_R1=r1, raw_R2=r2,
                                            barcode=barcode,
                                            sample_id=sample_id,
                                            threshold=threshold)
            log_file.write(command)

            process = Popen(command, stdout=PIPE, stderr=PIPE, shell=True)
            
            stdout, stderr = process.communicate()
            stdout_message = stdout.decode('ascii')
            stderr_message = stderr.decode('ascii')
            log_file.write(stdout_message)
            log_file.write(stderr_message)
            
            samples_done += 1
            
            if procedure == "mif2fastq":
                list_of_files.append((sample_id, "{}_R1.fastq".format(sample_id), "{}_R1.fastq".format(sample_id)))
            print_progress_bar(samples_done, samples_total, program_name=program_name)
    
    print("Written standard output to: {}".format(log_filename))
    if procedure == "mif2fastq":
        file_list_filename = os.path.join(output_folder, "assemble.log.txt")
        pd.DataFrame(list_of_files, columns=["#SAMPLE_ID", "OUTPUT_ASSEMBLY1", "OUTPUT_ASSEMBLY2"]).to_csv(file_list_filename, sep="\t", index=False)
        print("Written fastq file list to: {}".format(file_list_filename))  
    
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
    
def print_progress_bar(samples_done, samples_total, program_name="", objects="samples"):
    done = int(samples_done/samples_total*50)
    bar = '#'*done + '-'*(50-done)
    progress_bar = "{} |{}| {}/{} {} processed".format(program_name, bar, samples_done, samples_total, objects)
    end = "\r"
    if samples_done == samples_total:
        end = "\n"
    print(progress_bar, end=end)