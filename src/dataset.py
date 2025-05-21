import pandas as pd
import numpy as np
import torch
from torch.utils.data import Dataset


class DNASequenceDataset(Dataset):
    def __init__(self, csv_file):
        data = pd.read_csv(csv_file)
        self.sequences = [self.one_hot_encode(seq) for seq in data['Sequences']]
        self.labels = data['Labels'].values

    def __len__(self):
        return len(self.sequences)

    def __getitem__(self, idx):
        sequence = torch.tensor(self.sequences[idx], dtype=torch.float32)
        return sequence, torch.tensor(self.labels[idx], dtype=torch.long)

    def one_hot_encode(self, seq):
        mapping = {'A': [1, 0, 0, 0], 'C': [0, 1, 0, 0], 'G': [0, 0, 1, 0], 'T': [0, 0, 0, 1], 'N': [0, 0, 0, 0]}
        return np.array([mapping.get(base, [0, 0, 0, 0]) for base in seq])
