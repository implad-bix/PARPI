import os
import torch
import numpy as np
from src.model import CNN_Attention
import psutil
from datetime import datetime
import logging
from Bio import SeqIO
import gzip, gc
from torch.quantization import quantize_dynamic


def extended_fasta(input_fasta, output_fasta):
    with open(input_fasta, "r") as infile, open(output_fasta, "w") as outfile:
        for record in SeqIO.parse(infile, "fasta"):
            seq = str(record.seq)
            prefix = seq[:1000] if len(seq) >= 1000 else seq
            new_seq = seq + prefix
            new_record = record
            new_record.seq = new_seq
            SeqIO.write(new_record, outfile, "fasta")
    return output_fasta


def gfa_to_fasta(gfa_file, fasta_file):
    """
    Convert sequences from a GFA file to a FASTA file.

    Args:
    gfa_file (str): The path to the input GFA file.
    fasta_file (str): The path to the output FASTA file.
    """
    with open(gfa_file, 'r') as gfa:
        with open(fasta_file, 'w') as fasta:
            for line in gfa:
                # Skip header lines or empty lines
                if line.startswith('H') or not line.strip():
                    continue
                # Process the 'S' lines that contain sequences
                if line.startswith('S'):
                    # Split the line by tabs
                    parts = line.strip().split('\t')
                    # The first part is the node name (ID)
                    node_id = parts[1]
                    # The second part is the sequence
                    sequence = parts[2]
                    # Write the sequence to the FASTA file
                    fasta.write(f">{node_id}\n{sequence}\n")


def count_sequences(file_path, file_type):
    seq_count = 0
    with open(file_path, 'r') as file:
        lines = file.readlines()
    if file_type == "fastq":
        seq_count = sum(1 for line in lines if line.startswith('@'))
    elif file_type == "fasta":
        seq_count = sum(1 for line in lines if line.startswith('>'))
    return seq_count


def extract_reads(inputfile, datatype, data_size, outputfile):
    data_size = data_size * 1000
    total_length = 0
    sequences_to_save = []
    for record in SeqIO.parse(inputfile, datatype):
        sequence_length = len(record.seq)
        total_length += sequence_length
        sequences_to_save.append(record)
        if total_length > data_size:
            break
    with open(outputfile, "w") as output_handle:
        SeqIO.write(sequences_to_save, output_handle, datatype)
    return outputfile


def remove_duplicate_reads(input_file, file_format, output_dir='./'):
    log_file = os.path.join(output_dir, "ChloroScan.log.txt")
    logging.basicConfig(filename=log_file, level=logging.INFO, format='%(message)s', filemode='a', encoding='utf-8')
    current_date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    output_file = os.path.join(output_dir, f"plastid.rmdup.{file_format}")
    #remove duplicate reads ID#
    print(f"{current_date} - Removing duplicate reads.\n")
    logging.info(f"Removing duplicate reads.\n")
    seen_ids = set()
    with open(output_file, "w") as output_handle:
        for record in SeqIO.parse(input_file, file_format):
            if record.id in seen_ids:
                continue
            seen_ids.add(record.id)
            SeqIO.write(record, output_handle, file_format)
    current_date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"{current_date} - Removed duplicate reads based on reads ID and save it to {output_file}.")
    logging.info(f"Removed duplicate reads based on reads ID and save it to {output_file}.\n")

    return output_file


def one_hot_encode(seqs):
    """One-hot encode a list of DNA sequences"""
    mapping = {'A': [1, 0, 0, 0], 'C': [0, 1, 0, 0], 'G': [0, 0, 1, 0], 'T': [0, 0, 0, 1], 'N': [0, 0, 0, 0]}
    return np.array([np.array([mapping.get(base, [0, 0, 0, 0]) for base in seq]) for seq in seqs])


def get_memory_usage():
    """Get the current memory usage of the process (in MB)"""
    process = psutil.Process()
    return process.memory_info().rss / (1024 * 1024)


def load_model(model_path):
    """Load the trained model and apply quantization"""
    try:
        # Check if the model is a TorchScript model (.pt file)
        if model_path.endswith('.pt'):
            model = torch.jit.load(model_path)  # Load TorchScript model
        else:
            model = CNN_Attention(input_size=4, output_size=2)
            model.load_state_dict(torch.load(model_path, weights_only=True))  # Load weights for normal PyTorch model

            # Apply dynamic quantization to the model
            model = quantize_dynamic(model, {torch.nn.Linear, torch.nn.LSTM}, dtype=torch.qint8)  # Quantize Linear and LSTM layers
        
        model.eval()
        return model
    except Exception as e:
        print(f"Failed to load model: {e}")
        raise



def predict_batch(model, sequences):
    """Predict for a batch of sequences"""
    encoded_sequences = one_hot_encode(sequences)
    encoded_sequences = torch.tensor(encoded_sequences, dtype=torch.float32)

    with torch.no_grad():
        outputs = model(encoded_sequences)

    softmax = torch.nn.Softmax(dim=1)
    probs = softmax(outputs)
    return probs[:, 1].tolist()  # Return probability for class 1



def is_gzipped(file_path):
    """Check if a file is gzipped"""
    return file_path.endswith('.gz')


def open_file(file_path, datatype):
    """Open a file, handle both plain and gzip formats"""
    if is_gzipped(file_path):
        return gzip.open(file_path, 'rt', encoding='utf-8')
    else:
        return open(file_path, 'r', encoding='utf-8')


def sanitize_sequence(seq):
    """Sanitize a sequence by replacing invalid characters with 'N'"""
    return ''.join([base if base in 'ACGTN' else 'N' for base in seq])


def split_sequence(seq, segment_length=150, step_size=450):
    """Split a long sequence into smaller segments"""
    segments = [seq[i:i + segment_length] for i in range(0, len(seq) - segment_length + 1, step_size)]
    return segments


def is_gzip_file(file_path):
    with open(file_path, 'rb') as f:
        return f.read(2) == b'\x1f\x8b'


def process_single_sequences(data_size, input_file, datatype, target_length=1000, output_dir="./",
                                model_dir="./", batch_size=500, num_workers=4):
    if data_size == 'all':
        max_data_size = None
    else:
        max_data_size = int(data_size) * 1_000
    current_date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    output_file = os.path.join(output_dir, f"plastid.{datatype}")
    log_file = os.path.join(output_dir, "ChloroScan.log.txt")

    written_count = 0
    total_written_size = 0

    logging.basicConfig(filename=log_file, level=logging.INFO, format='%(message)s', filemode='a', encoding='utf-8')
    logging.info(f"Processing started on {current_date}\n")
    print(f"{current_date} - Processing started on {current_date}\n")

    model = load_model(model_dir)  # Load the model once at the beginning
    if is_gzip_file(input_file):
        open_fn = gzip.open
    else:
        open_fn = open
    with open_fn(input_file, 'rt') as handle, open(output_file, 'a') as output_handle:
        batch = []  # List to accumulate sequences for batch processing
        for i, record in enumerate(SeqIO.parse(handle, datatype)):
            seq = str(record.seq)
            seq = sanitize_sequence(seq)  # Sanitize the sequence

            # Process each sequence and add to the batch if the sequence length is sufficient
            subsequences = split_sequence(seq, segment_length=target_length, step_size=target_length)
            batch.append((record, subsequences))
            # If batch is full, process the batch
            if len(batch) == batch_size:
                subsequences_batch = [subseq for _, subseqs in batch for subseq in subseqs]
                #print('batch: '+str(i))
                # Predict probabilities for all subsequences in the batch
                probs = predict_batch(model, subsequences_batch)
                #print(len(subsequences_batch))
                #print('predict_batch: ' + str(i))
                # Process each parent sequence and calculate the weighted average probability
                for idx, (record, subsequences) in enumerate(batch):
                    subseq_probs = probs[:len(subsequences)]
                    probs = probs[len(subsequences):]

                    # Calculate the weighted average probability for this parent sequence
                    weighted_avg_prob = np.mean(subseq_probs)

                    # Check if the weighted average probability exceeds the threshold
                    if weighted_avg_prob > 0.5:
                        SeqIO.write(record, output_handle, datatype)
                        written_count += 1
                        total_written_size += len(record.seq)
                # Log progress after processing the batch
                if (i + 1) % 2000 == 0:
                    current_date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                    print_message = (f"{current_date} - Processed {i + 1} sequences. Written {written_count} "
                                    f"plastid sequences. Total written size: {total_written_size / 1_000_000:.2f} MB.")
                    progress_message = (f"Processed {i + 1} sequences. Written {written_count} "
                                        f"plastid sequences. Total written size: {total_written_size / 1_000_000:.2f} MB.")
                    print(print_message)
                    logging.info(progress_message)

                # Clear batch for next iteration
                batch.clear()

                # Stop if the total written size exceeds the limit
                if max_data_size != None:
                    if total_written_size >= max_data_size:
                        #print(f"Data size limit of {data_size / 1_000_000:.2f} MB reached. Stopping processing.")
                        break
                # Perform garbage collection to free up memory
                gc.collect()

        # If there are remaining sequences that don't make up a full batch, process them
        if batch:
            subsequences_batch = [subseq for _, subseqs in batch for subseq in subseqs]

            # Predict probabilities for the last batch
            probs = predict_batch(model, subsequences_batch)

            # Process each parent sequence and calculate the weighted average probability
            for idx, (record, subsequences) in enumerate(batch):
                subseq_probs = probs[:len(subsequences)]
                probs = probs[len(subsequences):]

                # Calculate the weighted average probability for this parent sequence
                weighted_avg_prob = np.mean(subseq_probs)

                # Check if the weighted average probability exceeds the threshold
                if weighted_avg_prob > 0.5:
                    SeqIO.write(record, output_handle, datatype)
                    written_count += 1
                    total_written_size += len(record.seq)

            # Log progress after processing the last batch
            current_date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            print_message = (f"{current_date} - Written {written_count} "
                                f"sequences. Total written size: {total_written_size / 1_000_000:.2f} MB.")
            progress_message = (f"Written {written_count} "
                                f"sequences. Total written size: {total_written_size / 1_000_000:.2f} MB.")
            print(print_message)
            logging.info(progress_message)

    # Final message once processing is complete
    current_date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print_final_message = (f"{current_date} - Selected plastid sequence count: {written_count}\n"
                     f"{current_date} - Total written size: {total_written_size / 1_000_000:.2f} MB.\n")
    final_message = (f"Selected plastid sequence count: {written_count}\t"
                     f"Total written size: {total_written_size / 1_000_000:.2f} MB.\n")
    logging.info(final_message)
    print(print_final_message)
    return output_file