import os


def count_sequences(file_path, file_type):
    seq_count = 0
    with open(file_path, 'r') as file:
        lines = file.readlines()
    if file_type == "fastq":
        seq_count = sum(1 for line in lines if line.startswith('@'))
    elif file_type == "fasta":
        seq_count = sum(1 for line in lines if line.startswith('>'))
    return seq_count


def process_fastq(reference_path, file_path, file_type, output_path, minimap2_path, seqkit_path, threads):
    paf_path = os.path.join(output_path, 'mapped.paf')
    select_paf_path = os.path.join(output_path, 'mapped_select.paf')
    readsID_path = os.path.join(output_path, 'mapped.readsID.txt')
    mapped_path = os.path.join(output_path, f'mapped.{file_type}')
    mapped_rmdup_path = os.path.join(output_path, f'mapped.rmdup.{file_type}')
    os.system(f'{minimap2_path} -c -t {threads} {reference_path} {file_path} > {paf_path}')
    count = 0
    with open(paf_path, "r") as file, open(select_paf_path, 'a+') as file2:
        for line in file:
            columns = line.strip().split("\t")
            if len(columns) > 0:
                reads_length = int(columns[1])
                reads_aligned = int(columns[9])
                ref_aligned = int(columns[10])
                if reads_aligned/reads_length >= 0.7 and reads_aligned/ref_aligned >= 0.7 and reads_length >= 1000:
                    file2.write(line)
                    count += 1
    data = []
    seen = set()
    with open(select_paf_path, "r") as file:
        for line in file:
            columns = line.strip().split("\t")
            if len(columns) > 0:
                first_column = columns[0]
                if first_column not in seen:
                    data.append(first_column)
                    seen.add(first_column)
    with open(readsID_path, "a+") as file:
        for item in data:
            file.write(item + '\n')
    os.system(f'{seqkit_path} grep -f {readsID_path} {file_path} > {mapped_path}')
    os.system(f'{seqkit_path} rmdup -n -i {mapped_path} -o {mapped_rmdup_path}')
    os.system(f'rm {mapped_path}')
    return mapped_rmdup_path, count

