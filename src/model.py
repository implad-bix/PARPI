import torch
import torch.nn as nn
import torch.nn.functional as F
import math


class PositionalEncoding(nn.Module):
    def __init__(self, d_model, max_len=5000):
        super().__init__()
        position = torch.arange(max_len).unsqueeze(1)
        div_term = torch.exp(torch.arange(0, d_model, 2) * (-math.log(10000.0) / d_model))
        pe = torch.zeros(1, max_len, d_model)
        pe[0, :, 0::2] = torch.sin(position * div_term)
        pe[0, :, 1::2] = torch.cos(position * div_term)
        self.register_buffer('pe', pe)

    def forward(self, x):
        x = x + self.pe[:, :x.size(1)]
        return x


class MultiScaleConv(nn.Module):
    def __init__(self, input_size):
        super().__init__()
        self.conv3 = nn.Conv1d(input_size, 32, kernel_size=3, padding=1)
        self.conv7 = nn.Conv1d(input_size, 32, kernel_size=7, padding=3)
        self.conv15 = nn.Conv1d(input_size, 32, kernel_size=15, padding=7)

    def forward(self, x):
        x3 = F.relu(self.conv3(x))
        x7 = F.relu(self.conv7(x))
        x15 = F.relu(self.conv15(x))
        return torch.cat([x3, x7, x15], dim=1)


class ResidualConvBlock(nn.Module):
    def __init__(self, in_channels):
        super().__init__()
        self.conv1 = nn.Conv1d(in_channels, in_channels * 2, kernel_size=5, padding=2)
        self.conv2 = nn.Conv1d(in_channels * 2, in_channels, kernel_size=3, padding=1)
        self.bn = nn.BatchNorm1d(in_channels)

    def forward(self, x):
        residual = x
        x = F.relu(self.conv1(x))
        x = self.conv2(x)
        x = self.bn(x + residual)
        return F.relu(x)


class MultiHeadAttention(nn.Module):
    def __init__(self, hidden_size, num_heads=8):
        super().__init__()
        self.num_heads = num_heads
        self.head_dim = hidden_size // num_heads

        self.query = nn.Linear(hidden_size, hidden_size)
        self.key = nn.Linear(hidden_size, hidden_size)
        self.value = nn.Linear(hidden_size, hidden_size)
        self.fc_out = nn.Linear(hidden_size, hidden_size)

    def forward(self, x):
        batch_size = x.shape[0]

        Q = self.query(x)
        K = self.key(x)
        V = self.value(x)

        Q = Q.view(batch_size, -1, self.num_heads, self.head_dim).permute(0, 2, 1, 3)
        K = K.view(batch_size, -1, self.num_heads, self.head_dim).permute(0, 2, 1, 3)
        V = V.view(batch_size, -1, self.num_heads, self.head_dim).permute(0, 2, 1, 3)

        energy = torch.matmul(Q, K.permute(0, 1, 3, 2)) / (self.head_dim ** 0.5)
        attention = torch.softmax(energy, dim=-1)

        x = torch.matmul(attention, V)
        x = x.permute(0, 2, 1, 3).contiguous()
        x = x.view(batch_size, -1, self.num_heads * self.head_dim)
        x = self.fc_out(x)
        return x


class CNN_Attention(nn.Module):
    def __init__(self, input_size, output_size, num_heads=8, dropout=0.3):
        super().__init__()

        self.pos_encoder = PositionalEncoding(d_model=4)

        self.conv_layers = nn.Sequential(
            MultiScaleConv(input_size),
            nn.MaxPool1d(2),
            ResidualConvBlock(96),
            nn.Conv1d(96, 256, kernel_size=5, padding=2),
            nn.ReLU(),
            nn.MaxPool1d(2),
            nn.Conv1d(256, 512, kernel_size=3, padding=1),
            nn.ReLU()
        )

        self.attn = MultiHeadAttention(hidden_size=512, num_heads=num_heads)

        self.classifier = nn.Sequential(
            nn.Linear(512, 256),
            nn.BatchNorm1d(256),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(256, 128),
            nn.ReLU(),
            nn.Linear(128, output_size)
        )

    def forward(self, x):
        x = self.pos_encoder(x)
        x = x.permute(0, 2, 1)

        x = self.conv_layers(x)
        x = x.permute(0, 2, 1)

        attn_out = self.attn(x)
        x = attn_out.mean(dim=1)

        return self.classifier(x)