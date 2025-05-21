import os, subprocess, sys, argparse, gzip, logging
from datetime import datetime
from src.processed_reads import extract_reads, count_sequences, gfa_to_fasta, remove_duplicate_reads, \
     process_single_sequences, extended_fasta
from src.minimap2mapping import process_fastq


def check_fasta(file_path):
    try:
        with open(file_path, 'r') as file:
            lines = file.readlines()
        if not lines:
            raise ValueError("The file is empty")
        for line in lines:
            if line.strip():
                if line.startswith(">"):
                    continue
                elif not all(char in 'ACGTacgt' for char in line.strip()):
                    raise ValueError(f"File content format error, not FASTA format: {line.strip()}")
            else:
                continue
        #print("The file is a valid FASTA format")
    except Exception as e:
        print(f"Error: {e}")


def check_file_format(file_path, data_type):
    """Check if the file format is correct (FASTA or FASTQ)"""
    is_gzipped = file_path.endswith(".gz")

    try:
        if is_gzipped:
            with gzip.open(file_path, 'rt') as f:
                content = f.readline()
        else:
            with open(file_path, 'r') as f:
                content = f.readline()
    except Exception as e:
        print(f"Error: Could not read the file {file_path}. {str(e)}")
        raise ValueError(f"Error: Could not read the file {file_path}. {str(e)}")

    if data_type == "fasta":
        if not content.startswith(">"):
            print(f"Error: {file_path} is not a valid FASTA file.")
            raise ValueError(f"Error: {file_path} is not a valid FASTA file.")
    elif data_type == "fastq":
        if not content.startswith("@"):
            print(f"Error: {file_path} is not a valid FASTQ file.")
            raise ValueError(f"Error: {file_path} is not a valid FASTQ file.")
    else:
        print("Error: Unsupported data type specified.")
        raise ValueError("Error: Unsupported data type specified.")


def validate_files(files, data_type):
    """Validate the input files based on data type"""
    for file_path in files:
        if not os.path.exists(file_path):
            print(f"Error: File {file_path} does not exist.")
            raise ValueError(f"Error: File {file_path} does not exist.")

        try:
            check_file_format(file_path, data_type)
        except ValueError as e:
            logging.error(e)
            sys.exit(1)


def calculate_data_size(depth, genome_size):
    """Calculate data size based on depth and genome size (in KB)"""
    if depth != 'all' and int(depth) < 50:
        print("Error: Depth must be at least 50 if a specific value is provided!")
        raise ValueError("Error: Depth must be at least 50 if a specific value is provided!")

    if depth == 'all':
        return 'all'  # Use all data without calculating the size
    else:
        # Calculate data size based on genome size (in KB) and depth
        return (genome_size) * int(depth)


def print_version():
    """Print the version number"""
    print("PARPI version 1.1")


def print_help():
    """Print help message"""
    help_text = """
    Usage:
    python PARPI.py --mode <1|2|3> --datatype <fasta|fastq> --file1 <file_path1> --file2 <file_path2> 
    --long_reads <file_path3> --genome_size <(int) kb> --depth <depth|all> --out_dir <output_directory> 
    --threads <num_threads> 
    [--help] [--version] [--reference] [--long_reads_type <type>]

    Arguments:
    --mode             Mode 1: Assembling long reads without error correction.
                       Mode 2: Assembling long reads, then short reads are used for error correction. 
                       Mode 3: Provide an assembled genome, then short reads are used for error correction.
                       (required)
    --datatype         Data type: either fasta or fastq (required)
    --long_reads       Path to the long reads file (required for mode 1/2)
    --long_reads_type  Type of long reads setting for flye: <pacbio-raw | pacbio-corr | pacbio-hifi | nano-raw |  
                       nano-corr | nano-hq> (required for mode 1/2, default: nano-raw)
    --file1            Path to the first file of paired-end reads (required for mode 2/3)
    --file2            Path to the second file of paired-end reads (required for mode 2/3)
    --assembled_genome Path to the assembled plastome (required for mode 3)
    --out_dir          Output directory path (required)
    --reference        Reference plastome file path (optional, default is 'no-reference')
    --genome_size      Plastome size in Kb (optional, default: 150, recommended range: 20-300)
    --depth            Coverage depth, must be an integer greater than or equal to 50, or 'all' (optional, default: 100)
    --threads          Number of threads to use for processing (optional, default: 10)

    --help             Display this help message
    --version          Display the version number

    Example:
    Mode 1:
    python PARPI.py --mode 1 --datatype fastq --long_reads long_reads.fastq --long_reads_type nano-raw 
    --genome_size 150 --depth 100 --reference path_to_ref_palstome.fasta --out_dir ./output --threads 10
    Mode 2:
    python PARPI.py --mode 2 --datatype fasta --long_reads long_reads.fasta --file1 reads_R1.fastq --file2 
    reads_R2.fastq --long_reads_type pacbio-hifi --round 3 --reference path_to_ref_palstome.fasta --genome_size 150 
    --depth all --out_dir ./output --threads 10
    Mode 3:
    python PARPI.py --mode 3 --datatype fastq --assembled_genome assembled_plastome.fasta --file1 reads_R1.fastq --file2 
    reads_R2.fastq --out_dir ./output --threads 10
    """
    print(help_text)


def validate_long_reads_type(long_reads_type):
    """Validate the long reads type"""
    valid_types = ["pacbio-raw", "pacbio-corr", "pacbio-hifi", "nano-raw", "nano-corr", "nano-hq"]
    if long_reads_type not in valid_types:
        print(f"Error: Invalid long_reads_type '{long_reads_type}'. Must be one of: pacbio-raw, pacbio-corr, pacbio-hifi, nano-raw, nano-corr, nano-hq.")
        sys.exit(1)


def using_flye(flye_path, long_reads_type, reads_send2flye, out_dir, threads, get_slim_graph_path):
    ################################################
    # Step 2: Using flye assembler to assemble plastome
    current_date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"{current_date} - #########################################################")
    print(f"{current_date} - ### Step 2: Using flye assembler to assemble plastome ###")
    print(f"{current_date} - #########################################################\n")
    logging.info(f"#########################################################")
    logging.info(f"### Step 2: Using flye assembler to assemble plastome ###")
    logging.info(f"#########################################################\n")
    print(f"{current_date} - Command line: {flye_path} --{long_reads_type} {reads_send2flye} "
          f"--out-dir {out_dir}/flye_plastome --threads {threads}\n")
    logging.info(f"Command line: {flye_path} --{long_reads_type} {reads_send2flye} "
                 f"--out-dir {out_dir}/flye_plastome --threads {threads}\n")
    os.system(
        f'{flye_path} --{long_reads_type} {reads_send2flye} --out-dir {out_dir}/flye_plastome --threads {threads}')
    current_date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    flye_gfa_dir = os.path.join(out_dir, 'flye_plastome', 'assembly_graph.gfa')
    # Check if the Flye result file exists
    if not os.path.exists(flye_gfa_dir):
        print(f"{current_date} - Flye assembler encountered an error, expected result not produced.\n")
        logging.info(f"Flye assembler encountered an error, expected result not produced.\n")
        raise FileNotFoundError("Flye assembler encountered an error, expected result not produced.")
    else:
        flye_gfa2fasta = os.path.join(out_dir, 'flye_plastome', 'assembly_graph.fasta')
        gfa_to_fasta(flye_gfa_dir, flye_gfa2fasta)
    current_date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"{current_date} - #######################")
    print(f"{current_date} - ### Step 2 finished ###")
    print(f"{current_date} - #######################\n")
    logging.info(f"#######################")
    logging.info(f"### Step 2 finished ###")
    logging.info(f"#######################\n")
    #####################################################
    # Step 3: Solve path of plastome from flye assembly results
    current_date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"{current_date} - #################################################################")
    print(f"{current_date} - ### Step 3: Solve path of plastome from flye assembly results ###")
    print(f"{current_date} - #################################################################\n")
    logging.info(f"#################################################################")
    logging.info(f"### Step 3: Solve path of plastome from flye assembly results ###")
    logging.info(f"#################################################################\n")
    flye_final_dir = os.path.join(out_dir, 'flye_final_assembly')
    print(f"{current_date} - Command line: {get_slim_graph_path} -F embplant_pt "
          f"-g {flye_gfa_dir} -o {flye_final_dir} -t {threads}\n")
    logging.info(f"Command line: {get_slim_graph_path} -F embplant_pt -g {flye_gfa_dir} "
                 f"-o {flye_final_dir} -t {threads}\n")
    os.system(f'{get_slim_graph_path} -F embplant_pt -g {flye_gfa_dir} -o {flye_final_dir} -t {threads}')
    # Check if any .fasta file exists in the directory
    fasta_files = [f for f in os.listdir(flye_final_dir) if f.endswith('.fasta')]
    current_date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    if not fasta_files:
        print(f"{current_date} - get_organelle_from_assembly.py encountered an error, expected result not produced.")
        logging.info(f"get_organelle_from_assembly.py encountered an error, expected result not produced.\n")
        sys.exit(1)
    else:
        long_reads_assembled_plastome = fasta_files[0]
        long_reads_assembled_plastome_path = os.path.join(flye_final_dir, long_reads_assembled_plastome)
    print(f"{current_date} - #######################")
    print(f"{current_date} - ### Step 3 finished ###")
    print(f"{current_date} - #######################\n")
    logging.info(f"#######################")
    logging.info(f"### Step 3 finished ###")
    logging.info(f"#######################\n")
    return long_reads_assembled_plastome_path, flye_gfa2fasta, flye_final_dir, fasta_files


def error_correction(out_dir, files, ropebwt2_path, msbwt_path, fasta_files, fmlrc2_path, flye_final_dir,
                     threads, round, extended_plastome, samtools_path, bwa_path, bedtools_path, datatype):
    # step 4 for mode 2 error correction
    ##################################################################
    # Step 4: Capture plastid reads from short reads data using minimap2
    current_date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    logging.info(f"##################################################")
    logging.info(f"### Step 4: Mapping short reads data using bwa ###")
    logging.info(f"##################################################\n")
    print(f"{current_date} - ##################################################")
    print(f"{current_date} - ### Step 4: Mapping short reads data using bwa ###")
    print(f"{current_date} - ##################################################\n")
    aligned_sam_dir = os.path.join(out_dir, 'aligned.sam')
    aligned_bam_dir = os.path.join(out_dir, 'aligned.bam')
    aligned_file1_dir = os.path.join(out_dir, f'plastome_filter_R1.{datatype}')
    aligned_file2_dir = os.path.join(out_dir, f'plastome_filter_R2.{datatype}')

    # bwa mapping short reads to ref_plastome
    print(f"{current_date} - {bwa_path} index {extended_plastome}\n")
    logging.info(f"Command line: {bwa_path} index {extended_plastome}\n")
    os.system(f'{bwa_path} index {extended_plastome}')
    print(
        f"{current_date} - Command line: {bwa_path} mem -t {threads} {extended_plastome} {files[1]} {files[2]} > {aligned_sam_dir}\n")
    logging.info(
        f"Command line: {bwa_path} mem -t {threads} {extended_plastome} {files[1]} {files[2]} > {aligned_sam_dir}\n")
    # os.system(f'{bwa_path} mem -t {threads} {long_reads_assembled_plastome_path} {files[1]} {files[2]} >
    # {aligned_sam_dir}')
    subprocess.run(f'{bwa_path} mem -t {threads} {extended_plastome} {files[1]} {files[2]} > {aligned_sam_dir}',
                   shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    os.system(f'rm {extended_plastome}*')
    # using samtools convert sam to bam
    current_date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"{current_date} - Command line: {samtools_path} view -@ {threads} -bS -F 12 {aligned_sam_dir} "
          f" > {aligned_bam_dir}\n")
    logging.info(f"Command line: {samtools_path} view -@ {threads} -bS -F 12 {aligned_sam_dir} "
                 f" > {aligned_bam_dir}\n")
    os.system(f'{samtools_path} view -@ {threads} -bS -F 12 {aligned_sam_dir} > {aligned_bam_dir}')
    os.system(f'rm {aligned_sam_dir}')
    # using bedtools convert bam to fastq/fasta
    current_date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"{current_date} - Command line: {bedtools_path} bamtofastq -i {aligned_bam_dir} -fq {aligned_file1_dir} "
          f" -fq2 {aligned_file2_dir}\n")
    logging.info(f"Command line: {bedtools_path} bamtofastq -i {aligned_bam_dir} -fq {aligned_file1_dir} "
                 f" -fq2 {aligned_file2_dir}\n")
    os.system(f'{bedtools_path} bamtofastq -i {aligned_bam_dir} -fq {aligned_file1_dir} -fq2 {aligned_file2_dir}')
    os.system(f'rm {aligned_bam_dir}')
    current_date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    if not os.path.exists(aligned_file1_dir) or not os.path.exists(aligned_file2_dir):
        logging.info(f"{current_date} - PARPI encountered an error, expected filtered plastid short reads "
                     f"not produced.\n")
        raise FileNotFoundError(
            "PARPI encountered an error, expected filtered plastid short reads not produced.")
    else:
        print(f"{current_date} - #######################")
        print(f"{current_date} - ### Step 4 finished ###")
        print(f"{current_date} - #######################\n")
        logging.info(f"#######################")
        logging.info(f"### Step 4 finished ###")
        logging.info(f"#######################\n")
        ##################################################################

    current_date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"{current_date} - ################################################################")
    print(f"{current_date} - ### Step 5: Error correction for long reads assembly results ###")
    print(f"{current_date} - ################################################################\n")
    logging.info(f"#####################################################################")
    logging.info(f"### Step Step 5: Error correction for long reads assembly results ###")
    logging.info(f"#####################################################################\n")
    polish_out_dir = os.path.join(out_dir, 'Polished_plastome')
    def is_compressed(file_path):
        return file_path.endswith(".gz")
    # Check if the files are compressed and execute the appropriate command
    if is_compressed(aligned_file1_dir):
        command = f"gunzip -c {aligned_file1_dir} {aligned_file2_dir}  | awk 'NR % 4 == 2' | sort " \
                    f"| tr NT TN | {ropebwt2_path} -LR | tr NT TN | {msbwt_path} convert {polish_out_dir}"
    else:
        command = f"cat {aligned_file1_dir} {aligned_file2_dir} | awk 'NR % 4 == 2' | sort | tr NT TN | {ropebwt2_path}  -LR | " \
                    f"tr NT TN | {msbwt_path} convert {polish_out_dir}"
    # Execute the command
    subprocess.run(command, shell=True, check=True)

    file_num = 1
    for file in fasta_files:
        input_file = os.path.join(flye_final_dir, file)
        for i in range(round):
            output_file = os.path.join(polish_out_dir, f'Polished{i + 1}_plastome_path{file_num}.fasta')

            subprocess.run(
                f'{fmlrc2_path} -t {threads} {polish_out_dir}/comp_msbwt.npy {input_file} '
                f'{output_file}',
                shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            input_file = output_file
        file_num += 1

    os.system(f'rm {polish_out_dir}/comp_msbwt.npy')
    current_date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"{current_date} - #######################")
    print(f"{current_date} - ### Step 5 finished ###")
    print(f"{current_date} - #######################\n")
    logging.info(f"#######################")
    logging.info(f"### Step 5 finished ###")
    logging.info(f"#######################\n")
    current_date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    logging.info(f"The final assembled plastome can be found in {polish_out_dir}\n")
    print(f"{current_date} - The final assembled plastome can be found in {polish_out_dir}\n")
    print(f"{current_date} - ######################################")
    print(f"{current_date} - ### PARPI complete assembly ###")
    print(f"{current_date} - ######################################\n")
    logging.info(f"######################################")
    logging.info(f"### PARPI complete assembly ###")
    logging.info(f"######################################\n")


def only_correction(out_dir, files, ropebwt2_path, msbwt_path, fasta_file, fmlrc2_path, threads, round, samtools_path,
                    bwa_path, bedtools_path, datatype):
    ##################################################################
    # Step 1: Mapping short reads data using bwa
    current_date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    logging.info(f"##################################################")
    logging.info(f"### Step 1: Mapping short reads data using bwa ###")
    logging.info(f"##################################################\n")
    print(f"{current_date} - ##################################################")
    print(f"{current_date} - ### Step 1: Mapping short reads data using bwa ###")
    print(f"{current_date} - ##################################################\n")
    aligned_sam_dir = os.path.join(out_dir, 'aligned.sam')
    aligned_bam_dir = os.path.join(out_dir, 'aligned.bam')
    aligned_file1_dir = os.path.join(out_dir, f'plastome_filter_R1.{datatype}')
    aligned_file2_dir = os.path.join(out_dir, f'plastome_filter_R2.{datatype}')

    # bwa mapping short reads to ref_plastome
    print(f"{current_date} - {bwa_path} index {fasta_file}\n")
    logging.info(f"Command line: {bwa_path} index {fasta_file}\n")
    os.system(f'{bwa_path} index {fasta_file}')
    print(
        f"{current_date} - Command line: {bwa_path} mem -t {threads} {fasta_file} {files[0]} {files[1]} > "
        f"{aligned_sam_dir}\n")
    logging.info(
        f"Command line: {bwa_path} mem -t {threads} {fasta_file} {files[0]} {files[1]} > {aligned_sam_dir}\n")
    # os.system(f'{bwa_path} mem -t {threads} {long_reads_assembled_plastome_path} {files[1]} {files[2]} >
    # {aligned_sam_dir}')
    subprocess.run(f'{bwa_path} mem -t {threads} {fasta_file} {files[0]} {files[1]} > {aligned_sam_dir}',
                   shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    os.system(f'rm {fasta_file}.*')
    # using samtools convert sam to bam
    current_date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"{current_date} - Command line: {samtools_path} view -@ {threads} -bS -F 12 {aligned_sam_dir} "
          f" > {aligned_bam_dir}\n")
    logging.info(f"Command line: {samtools_path} view -@ {threads} -bS -F 12 {aligned_sam_dir} "
                 f" > {aligned_bam_dir}\n")
    os.system(f'{samtools_path} view -@ {threads} -bS -F 12 {aligned_sam_dir} > {aligned_bam_dir}')
    os.system(f'rm {aligned_sam_dir}')
    # using bedtools convert bam to fastq/fasta
    current_date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"{current_date} - Command line: {bedtools_path} bamtofastq -i {aligned_bam_dir} -fq {aligned_file1_dir} "
          f" -fq2 {aligned_file2_dir}\n")
    logging.info(f"Command line: {bedtools_path} bamtofastq -i {aligned_bam_dir} -fq {aligned_file1_dir} "
                 f" -fq2 {aligned_file2_dir}\n")
    os.system(f'{bedtools_path} bamtofastq -i {aligned_bam_dir} -fq {aligned_file1_dir} -fq2 {aligned_file2_dir}')
    os.system(f'rm {aligned_bam_dir}')
    current_date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    if not os.path.exists(aligned_file1_dir) or not os.path.exists(aligned_file2_dir):
        logging.info(f"{current_date} - PARPI encountered an error, expected filtered plastid short reads "
                     f"not produced.\n")
        raise FileNotFoundError(
            "PARPI encountered an error, expected filtered plastid short reads not produced.")
    else:
        print(f"{current_date} - #######################")
        print(f"{current_date} - ### Step 1 finished ###")
        print(f"{current_date} - #######################\n")
        logging.info(f"#######################")
        logging.info(f"### Step 1 finished ###")
        logging.info(f"#######################\n")

    current_date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"{current_date} - ####################################################")
    print(f"{current_date} - ### Step 2 Error correction for specified genome ###")
    print(f"{current_date} - ####################################################\n")
    logging.info(f"####################################################")
    logging.info(f"### Step 2 Error correction for specified genome ###")
    logging.info(f"####################################################\n")
    polish_out_dir = os.path.join(out_dir, 'Polished_plastome')
    def is_compressed(file_path):
        return file_path.endswith(".gz")
    # Check if the files are compressed and execute the appropriate command
    if is_compressed(aligned_file1_dir):
        command = f"gunzip -c {aligned_file1_dir} {aligned_file2_dir}  | awk 'NR % 4 == 2' | sort " \
                    f"| tr NT TN | {ropebwt2_path} -LR | tr NT TN | {msbwt_path} convert {polish_out_dir}"
    else:
        command = f"cat {aligned_file1_dir} {aligned_file2_dir} | awk 'NR % 4 == 2' | sort | tr NT TN | " \
                  f"{ropebwt2_path}  -LR | " \
                    f"tr NT TN | {msbwt_path} convert {polish_out_dir}"
    # Execute the command
    subprocess.run(command, shell=True, check=True)

    input_file = fasta_file
    for i in range(round):
        output_file = os.path.join(polish_out_dir, f'Polished{i + 1}_plastome_path.fasta')
        subprocess.run(
            f'{fmlrc2_path} -t {threads} {polish_out_dir}/comp_msbwt.npy {input_file} '
            f'{output_file}',
            shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        input_file = output_file

    os.system(f'rm {polish_out_dir}/comp_msbwt.npy')

    current_date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    logging.info(f"The final polished plastome can be found in {polish_out_dir}\n")
    print(f"{current_date} - The final polished plastome can be found in {polish_out_dir}\n")
    print(f"{current_date} - ##############################################")
    print(f"{current_date} - ### PARPI complete error correction ###")
    print(f"{current_date} - ##############################################\n")
    logging.info(f"##############################################")
    logging.info(f"### PARPI complete error correction ###")
    logging.info(f"##############################################\n")


def process_mode(files, mode, datatype, long_reads_type, data_size, reference, threads, out_dir, flye_path,
                 minimap2_path, seqkit_path, get_slim_graph_path, model_path, ropebwt2_path, fmlrc2_path, msbwt_path,
                 round, assembled_genome, samtools_path, bwa_path, bedtools_path):
    ##################################################################
    # Step 1: Capture plastid reads from long reads data using minimap2
    if mode != '3':
        if reference != "no-reference":
            current_date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            print(f"{current_date} - ######################################################")
            print(f"{current_date} - ### Step 1: Mapping long reads data using minimap2 ###")
            print(f"{current_date} - ######################################################\n")
            logging.info(f"######################################################")
            logging.info(f"### Step 1: Mapping long reads data using minimap2 ###")
            logging.info(f"######################################################\n")
            print(f"{current_date} - The reference plastome is: {reference}\n")
            logging.info(f"The reference plastome is: {reference}\n")
            mapped_rmdup_path, count = process_fastq(reference, files[0], datatype, out_dir, minimap2_path,
                                                     seqkit_path, threads)
            if data_size == 'all':
                reads_send2flye = mapped_rmdup_path
                current_date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                logging.info(f"Mapped {count} potential plastid reads, sampled {count} potential plastid reads\n")
                print(f"{current_date} - Mapped {count} potential plastid reads, sampled {count} potential plastid reads\n")
            else:
                sampled_file = os.path.join(out_dir, f'mapped.sample.{datatype}')
                reads_send2flye = extract_reads(mapped_rmdup_path, datatype, data_size, sampled_file)
                seq_count = count_sequences(sampled_file, datatype)
                current_date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                logging.info(f"Mapped {count} potential plastid reads, sampled {seq_count} potential plastid reads\n")
                print(f"{current_date} - Mapped {count} potential plastid reads, sampled {seq_count} potential plastid reads\n")
            current_date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            print(f"{current_date} - #######################")
            print(f"{current_date} - ### Step 1 finished ###")
            print(f"{current_date} - #######################\n")
            logging.info(f"#######################")
            logging.info(f"### Step 1 finished ###")
            logging.info(f"#######################\n")
            ###step 2 and 3###
            long_reads_assembled_plastome_path, flye_gfa2fasta, \
            flye_final_dir, fasta_files = using_flye(flye_path, long_reads_type, reads_send2flye, out_dir, threads,
                                                     get_slim_graph_path)
            if mode == '1':
                current_date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                logging.info(f"The final assembled plastome can be found in {flye_final_dir}\n")
                print(f"{current_date} - The final assembled plastome can be found in {flye_final_dir}\n")
                print(f"{current_date} - ######################################")
                print(f"{current_date} - ### PARPI complete assembly ###")
                print(f"{current_date} - ######################################\n")
                logging.info(f"######################################")
                logging.info(f"### PARPI complete assembly ###")
                logging.info(f"######################################\n")
            elif mode == '2':
                ###step 4-5 ###
                extended_plastome = os.path.join(out_dir, f'extended_ref.fasta')
                extended_plastome = extended_fasta(long_reads_assembled_plastome_path, extended_plastome)
                error_correction(out_dir, files, ropebwt2_path, msbwt_path, fasta_files, fmlrc2_path, flye_final_dir,
                                 threads, round, extended_plastome, samtools_path, bwa_path, bedtools_path, datatype)
        else:
            current_date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            print(f"{current_date} - ########################################################")
            print(f"{current_date} - ### Step 1: Capture plastid reads based on CNN model ###")
            print(f"{current_date} - ########################################################\n")
            logging.info(f"########################################################")
            logging.info(f"### Step 1: Capture plastid reads based on CNN model ###")
            logging.info(f"########################################################\n")
            plastid_raw_reads = process_single_sequences(data_size, files[0], datatype, target_length=1000,
                                                         output_dir=out_dir, model_dir=model_path, batch_size=2)
            reads_send2flye = remove_duplicate_reads(plastid_raw_reads, datatype, output_dir=out_dir)
            current_date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            print(f"{current_date} - #######################")
            print(f"{current_date} - ### Step 1 finished ###")
            print(f"{current_date} - #######################\n")
            logging.info(f"#######################")
            logging.info(f"### Step 1 finished ###")
            logging.info(f"#######################\n")
            ###step 2 and 3###
            long_reads_assembled_plastome_path, flye_gfa2fasta, \
            flye_final_dir, fasta_files = using_flye(flye_path, long_reads_type, reads_send2flye, out_dir, threads,
                                                     get_slim_graph_path)
            if mode == '1':
                current_date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                logging.info(f"The final assembled plastome can be found in {flye_final_dir}\n")
                print(f"{current_date} - The final assembled plastome can be found in {flye_final_dir}\n")
                print(f"{current_date} - ######################################")
                print(f"{current_date} - ### PARPI complete assembly ###")
                print(f"{current_date} - ######################################\n")
                logging.info(f"######################################")
                logging.info(f"### PARPI complete assembly ###")
                logging.info(f"######################################\n")
            elif mode == '2':
                ###step 4-5 ###
                extended_plastome = os.path.join(out_dir, f'extended_ref.fasta')
                extended_plastome = extended_fasta(long_reads_assembled_plastome_path, extended_plastome)
                error_correction(out_dir, files, ropebwt2_path, msbwt_path, fasta_files, fmlrc2_path, flye_final_dir,
                                 threads, round, extended_plastome, samtools_path, bwa_path, bedtools_path, datatype)
    else:
        if assembled_genome:
            only_correction(out_dir, files, ropebwt2_path, msbwt_path, assembled_genome, fmlrc2_path,
                            threads, round, samtools_path, bwa_path, bedtools_path, datatype)


def main():
    # Pre-argument parsing to handle help/version flags
    if '--help' in sys.argv:
        print_help()
        return  # Return after printing help, do not continue further

    if '--version' in sys.argv:
        print_version()
        return  # Return after printing version, do not continue further

    # Setup argument parser
    parser = argparse.ArgumentParser()

    # Define arguments
    parser.add_argument("--mode", required=True, choices=["1", "2", "3"],
                        help="Mode 1: Assembling long reads without polish. Mode 2: "
                             "Assembling long reads, then short reads are used for error correction. "
                             "Mode 3: Provide an assembled genome, then short reads are used for error correction. "
                             "(required)")
    parser.add_argument("--datatype", required=True, choices=["fasta", "fastq"],
                        help="Input data type, either fasta or fastq (required)")
    parser.add_argument("--long_reads", help="Path to the long reads file (required for mode 1/2)")
    parser.add_argument("--long_reads_type", choices=["pacbio-raw", "pacbio-corr", "pacbio-hifi", "nano-raw", "nano-corr", "nano-hq"],
                        help="Type of long reads (required for mode 1/2)", default='nano-raw')
    parser.add_argument("--file1", help="Path to the first file of paired-end reads (required for mode 2/3)")
    parser.add_argument("--file2", help="Path to the second file of paired-end reads (required for mode 2/3)")
    #parser.add_argument("--round", type=int, default=3, help="Rounds for error correction (required for mode 2/3,
    # default: 3)")
    parser.add_argument("--assembled_genome", default=None, help="Path to the assembled plastome (required for mode 3)")
    parser.add_argument("--out_dir", required=True, help="Output directory path (required)")
    parser.add_argument("--reference", default="no-reference", help="Reference plastome file path (optional, default: 'no-reference')")
    parser.add_argument("--genome_size", type=int, default=150, help="Plastome size in Kb (optional, default: 150, recommended range: 20-300)")
    parser.add_argument("--depth", default="50", help="Coverage depth, must be an integer >= 50, or 'all' (optional, default: 50)")
    parser.add_argument("--threads", type=int, default=10, help="Number of threads to use for processing (optional, default: 10)")
    args = parser.parse_args()

    '''if args.round < 1:
        print("Error: --round must be an integer >= 1.")
        sys.exit(1)'''

    files = []
    # Validate input files based on mode
    if args.mode == "1":
        # Mode 1 requires both file1 and file2
        if args.file1 or args.file2 or not args.long_reads:
            print("Error: Mode 1 only require long reads data, and paired-end reads (file1, file2) are not required.")
            sys.exit(1)
        files = [args.long_reads]
        if args.assembled_genome:
            print("Error: --assembled_genome are not required for mode 1.")
            sys.exit(1)

    elif args.mode == "2":
        # Mode 2 requires only long_reads, file1 and file2 should not be present
        if not args.file1 or not args.file2 or not args.long_reads:
            print("Error: Mode 2 require long reads for assembling. While paired-end reads (file1, file2) "
                  "should be provided for error correction.")
            sys.exit(1)
        files = [args.long_reads, args.file1, args.file2]
        if args.assembled_genome:
            print("Error: --assembled_genome are not required for mode 2.")
            sys.exit(1)

    elif args.mode == "3":
        # Mode 2 requires only long_reads, file1 and file2 should not be present
        if not args.file1 or not args.file2 or args.long_reads:
            print("Error: Mode 3 require paired-end reads (file1, file2) for error correction. While long reads data ("
                  "file3) should not be provided.")
            sys.exit(1)
        files = [args.file1, args.file2]
        if not args.assembled_genome:
            print("Error: --assembled_genome are required for mode 3.")
            sys.exit(1)
        else:
            check_fasta(args.assembled_genome)
    # Validate files
    validate_files(files, args.datatype)

    # Calculate data size
    try:
        data_size = calculate_data_size(args.depth, args.genome_size)
    except ValueError as e:
        logging.error(e)
        sys.exit(1)

    # Handle output directory
    if args.out_dir:
        out_dir_abs_path = os.path.abspath(args.out_dir)
        if not os.path.exists(out_dir_abs_path):
            os.makedirs(out_dir_abs_path)
        logging.basicConfig(filename=os.path.join(out_dir_abs_path, "PARPI.log.txt"),
                            level=logging.INFO, format='%(asctime)s - %(message)s')
        current_date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print(f"{current_date} - Output directory: {out_dir_abs_path}\n")
        logging.info(f"Output directory: {out_dir_abs_path}\n")
    else:
        print("Error: You must specify an output directory (--out-dir).")
        sys.exit(1)

    logging.info(f"Starting PARPI with the following parameters: "
                 f"Mode: {args.mode}, Data type: {args.datatype}, "
                 f"Depth: {args.depth}, Genome size: {args.genome_size} KB, "
                 f"Output directory: {out_dir_abs_path}, Threads: {args.threads}, "
                 f"Reference: {args.reference}, Files: {files}\n")

    # Continue with main logic (Mode 1 or Mode 2 processing)
    current_script_path = os.path.dirname(os.path.realpath(__file__))
    get_slim_graph_path = os.path.join(current_script_path, 'src', 'get_organelle_from_assembly.py')
    minimap2_path = os.path.join(current_script_path, 'apps', 'minimap2')
    model_path = os.path.join(current_script_path, 'models', 'best_model.pth')
    #makeblastdb_path = os.path.join(current_script_path, 'apps', 'makeblastdb')
    #BLASTn_path = os.path.join(current_script_path, 'apps', 'blastn')
    seqkit_path = os.path.join(current_script_path, 'apps', 'seqkit')
    samtools_path = os.path.join(current_script_path, 'apps', 'samtools')
    bedtools_path = os.path.join(current_script_path, 'apps', 'bedtools')
    bwa_path = os.path.join(current_script_path, 'apps', 'bwa')
    ropebwt2_path = os.path.join(current_script_path, 'apps', 'ropebwt2')
    fmlrc2_path = os.path.join(current_script_path, 'apps', 'fmlrc2')
    flye_path = os.path.join(current_script_path, 'apps', 'flye-master', 'bin','flye')
    #spades_path = '/home/user001/mambaforge/envs/jll/bin//spades.py'
    msbwt_path = os.path.join(current_script_path, 'apps', 'msbwt-master', 'bin','msbwt')
    round = 1
    process_mode(files, args.mode, args.datatype, args.long_reads_type, data_size, args.reference, args.threads,
                 out_dir_abs_path, flye_path, minimap2_path, seqkit_path, get_slim_graph_path, model_path,
                 ropebwt2_path, fmlrc2_path, msbwt_path, round, args.assembled_genome,
                 samtools_path, bwa_path, bedtools_path)


if __name__ == "__main__":
    main()
