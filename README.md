# PARPI
You can download the compressed file ​​PARPI.zip​​ through the release link (https://github.com/implad-bix/PARPI/releases/tag/PARPI). 
After extracting it, ensure all files have the necessary permissions for use.

Before using PARPI, please ensure that PyTorch is installed, which can be done via Mamba:

mamba install pytorch -c pytorch.

Usage:
    python PARPI.py --mode <1|2|3> --datatype <fasta|fastq> --file1 <file_path1> --file2 <file_path2>
    --long_reads <file_path3> --genome_size <(int) kb> --depth <depth|all> --out_dir <output_directory>
    --threads <num_threads>
    [--help] [--version] [--reference] [--long_reads_type <type>]

    Arguments:
    --mode                      Mode 1: Assembling long reads without error correction.
                                Mode 2: Assembling long reads, then short reads are used for error correction.
                                Mode 3: Provide an assembled genome, then short reads are used for error correction.
                                (required)
    --datatype                  Data type: either fasta or fastq (required)
    --long_reads                Path to the long reads file (required for mode 1/2)
    --long_reads_type           Type of long reads setting for flye: <pacbio-raw | pacbio-corr | pacbio-hifi | nano-raw |
                                                                      nano-corr | nano-hq> (required for mode 1/2, default: nano-raw)
    --file1                     Path to the first file of paired-end reads (required for mode 2/3)
    --file2                     Path to the second file of paired-end reads (required for mode 2/3)
    --assembled_genome          Path to the assembled plastome (required for mode 3)
    --out_dir                   Output directory path (required)
    --reference                 Reference plastome file path (optional, default is 'no-reference')
    --genome_size               Plastome size in Kb (optional, default: 150, recommended range: 20-300)
    --depth                     Coverage depth, must be an integer greater than or equal to 50, or 'all' (optional, default: 100)
    --threads                   Number of threads to use for processing (optional, default: 10)

    --help                      Display this help message
    --version                   Display the version number

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
