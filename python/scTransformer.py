##Import the required packages
import torch
import torch.nn as nn
import warnings
warnings.filterwarnings('ignore')
import math
import pandas as pd
import scanpy as sc
import os
from tqdm.auto import tqdm
from torch.utils.data import (DataLoader,Dataset)
torch.set_default_tensor_type(torch.DoubleTensor)
import numpy as np
import random
from sklearn import preprocessing
##Set random seeds
def same_seeds(seed):
    random.seed(seed)
    # Numpy
    np.random.seed(seed)
    # Torch
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.benchmark = False
    torch.backends.cudnn.deterministic = True
same_seeds(2021)
class myDataSet(Dataset):
    def __init__(self, data, label):
        self.data = data
        self.label = label
        self.length = len(data)

    def __len__(self):
        return self.length

    def __getitem__(self, index):
        data = torch.from_numpy(self.data)
        label = torch.from_numpy(self.label)

        return data[index], label[index]
##Positional Encoder Layer
class PositionalEncoding(nn.Module):
    def __init__(self, d_model, dropout=0.1, max_len=5000):
        super(PositionalEncoding, self).__init__()
        self.dropout = nn.Dropout(p=dropout)
        pe = torch.zeros(max_len, d_model)
        position = torch.arange(0, max_len, dtype=torch.float).unsqueeze(1)
        div_term = torch.exp(torch.arange(0, d_model, 2).float() * (-math.log(10000.0) / d_model))
        
        ##the sine function is used to represent the odd-numbered sub-vectors
        pe[:, 0::2] = torch.sin(position * div_term)
        ##the cosine function is used to represent the even-numbered sub-vectors
        pe[:, 1::2] = torch.cos(position * div_term)
        pe = pe.unsqueeze(0).transpose(0, 1)
        self.register_buffer('pe', pe)

    def forward(self, x):
        x = x + self.pe[:x.size(0), :]
        return self.dropout(x)
class scTransformer(nn.Module):
    def __init__(self, input_dim, nhead=2, d_model=80, num_classes=2, dropout=0.1):
        super().__init__()
        #TransformerEncoderLayer with self-attention 
        self.encoder_layer = nn.TransformerEncoderLayer(
            d_model=d_model, dim_feedforward=1024, nhead=nhead, dropout=dropout
        )
        
        #Positional Encoding with self-attention 
        self.positionalEncoding = PositionalEncoding(d_model=d_model, dropout=dropout)
        
        #Classification layer
        self.pred_layer = nn.Sequential(
            nn.Linear(d_model, d_model),
            nn.ReLU(),
            nn.Linear(d_model, num_classes)
        )

    def forward(self, mels):
        out = mels.permute(1, 0, 2)
        #Positional Encoding layer
        out = self.positionalEncoding(out)
        #Transformer Encoder layer layer
        out = self.encoder_layer(out)
        out = out.transpose(0, 1)
        #Pooling layer
        out = out.mean(dim=1)
        #Classification layer
        out = self.pred_layer(out)
        return out
#Set parameters of scTransformer
d_models = 1024 
heads = 64                     #the number of heads in self-attention mechanism
num_classes = len(cell_types)  #the number of cell types

lr = 0.0001                    #learning rate
dp = 0.1                       #dropout rate
n_epochs = 50                  #the number of epoch
    
#Constructing the scTransformer model
model =  scTransformer(input_dim=d_models, nhead=heads, d_model=d_models,
                       num_classes=num_classes,dropout=dp)
#Setting loss function
criterion = nn.CrossEntropyLoss()
#Setting optimization function
optimizer = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=1e-5)
