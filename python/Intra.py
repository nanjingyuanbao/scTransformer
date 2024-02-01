##Import the required packages
import torch
import torch.nn as nn
import warnings
warnings.filterwarnings('ignore')
import os
from sklearn.metrics import accuracy_score,f1_score,recall_score,precision_score
import math
import pandas as pd
import scanpy as sc
from tqdm.auto import tqdm
from torch.utils.data import (DataLoader,Dataset)
torch.set_default_tensor_type(torch.DoubleTensor)
import numpy as np
import random
from sklearn.model_selection import KFold, StratifiedKFold
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
adata = sc.read_csv("/home/yuanjiaxin/ccc/xin/xinpancreas1.csv", delimiter=",")
if True:
        adata = adata.T
sc.pp.filter_genes(adata, min_cells=1)
#Preprocessing the scRNA-seq data
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
#Obtaining the HVGs
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
adata.raw = adata
adata = adata[:, adata.var.highly_variable]
# Converting the gene expression matrix into sub-vectors
#(n_cells,n_genes) -> (n_cells,gap_num,gap)  gap_num = int(gene_num / gap) + 1
single_cell_list = []
for single_cell in X:
        feature = []
        length = len(single_cell)
        
        gap=1024
        #spliting the gene expression vector into some sub-vectors whose length is gap
        for k in range(0, length, gap):
            if (k + gap > length):
                a = single_cell[length - gap:length]
            else:
                a = single_cell[k:k + gap]

            #scaling each sub-vectors 
            a = preprocessing.scale(a,axis=0, with_mean=True, with_std=True, copy=True)
            feature.append(a)

        feature = np.asarray(feature)
        single_cell_list.append(feature)
y_train = pd.read_csv("/home/yuanjiaxin/ccc/xin/Labels.csv")
y_train = y_train.T
y_train = y_train.values[0]
cell_types = []
labelss = []
for i in y_train:
    i = str(i).upper()
    if (not cell_types.__contains__(i)):
                cell_types.append(i)
    labelss.append(cell_types.index(i))
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
from sklearn.model_selection import KFold, StratifiedKFold
from sklearn.metrics import accuracy_score,f1_score,recall_score,precision_score
from tqdm.auto import tqdm
Data, labels, cell_types = single_cell_list, labelss, cell_types
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

l = Data.shape[0]
#Setting the batch size varies according to the number of cells
if(l > 5000):
        batch_sizes = 512   #the size of batch
else:
        batch_sizes = 256   #the size of batch
#Model verification by 5-KFold
skf = StratifiedKFold(n_splits=5,random_state=2021, shuffle=True)
fold = 0
Indexs = []
for index in range(len(Data)):
        Indexs.append(index)
Indexs = np.asarray(Indexs)
for train_index, test_index in skf.split(Data, labels):
        fold = fold + 1
        print(train_index)
        #print(test_index)
        X_train, X_test = Data[train_index], Data[test_index]
        X_train = np.asarray(X_train)
        X_test = np.asarray(X_test)
        
        y_train, y_test = labels[train_index], labels[test_index]
        
        y_train = np.asarray(y_train)
        y_test = np.asarray(y_test)
        print(X_train.shape)
        print(X_test.shape)
        print(y_train.shape)
        print(y_test.shape)
for train_index, test_index in skf.split(Data, labels):
        fold = fold + 1
        X_train, X_test = Data[train_index], Data[test_index]
        X_train = np.asarray(X_train)
        X_test = np.asarray(X_test)

        y_train, y_test = labels[train_index], labels[test_index]

        y_train = np.asarray(y_train)
        y_test = np.asarray(y_test)
        
        #Setting the training dataset
        train_dataset = myDataSet(data=X_train,label=y_train)
        train_loader = DataLoader(train_dataset,batch_size=batch_sizes,shuffle=True,pin_memory=True)
        
        #Setting the test dataset
        test_dataset = myDataSet(data=X_test,label=y_test)
        test_loader = DataLoader(test_dataset,batch_size=batch_sizes,shuffle=False,pin_memory=True)
        #startting training scTransformer.Using training data to train scTransformer
        #n_epochs: the times of Training 
        model.train()
        for epoch in range(n_epochs):
            # model.train()
            # These are used to record information in training.
            train_loss = []
            train_accs = []
            train_f1s = []
            for batch in tqdm(train_loader):
                # A batch consists of scRNA-seq data and corresponding cell type annotations.
                data, labels = batch
                logits = model(data)
                labels = torch.tensor(labels, dtype=torch.long)
                loss = criterion(logits, labels)
                optimizer.zero_grad()
                loss.backward()
                optimizer.step()
                #Getting the predicted cell type
                preds = logits.argmax(1)
                preds = preds.cpu().numpy()
                labels = labels.cpu().numpy()
                #Metrics
                acc = accuracy_score(labels, preds)
                f1 = f1_score(labels,preds,average='macro')
                train_loss.append(loss.item())
                train_accs.append(acc)
                train_f1s.append(f1)
            train_loss = sum(train_loss) / len(train_loss)
            train_acc = sum(train_accs) / len(train_accs)
            train_f1 = sum(train_f1s) / len(train_f1s)

            print(f"[ Train | {epoch + 1:03d}/{n_epochs:03d} ] loss = {train_loss:.5f}, acc = {train_acc:.5f}, f1 = {train_f1:.5f}")
            #print(f"[ Train | {epoch + 1:03d}/{n_epochs:03d} ] loss = {train_loss:.5f}, acc = {train_acc:.5f}, f1 = {train_f1:.5f}")
            ##Start the validation model, which predicts the cell types in the test dataset
            model.eval()
            test_accs = []
            test_f1s = []
            y_predict = []
            labelss = []
            for batch in tqdm(test_loader):
                # A batch consists of scRNA-seq data and corresponding cell type annotations.
                data, labels = batch
                with torch.no_grad():
                    logits = model(data)
                    
                #Getting the predicted cell type
                preds = logits.argmax(1)
                preds = preds.cpu().numpy().tolist()
                labels = labels.cpu().numpy().tolist()
                
                #Metrics
                acc = accuracy_score(labels, preds)
                f1 = f1_score(labels, preds, average='macro')
                test_f1s.append(f1)
                test_accs.append(acc)
                
                y_predict.extend(preds)
                labelss.extend(labels)
            test_acc = sum(test_accs) / len(test_accs)
            test_f1 = sum(test_f1s) / len(test_f1s)
            print("---------------------------------------------end test---------------------------------------------")
            #Metrics
            all_acc = accuracy_score(labelss, y_predict)
            all_f1 = f1_score(labelss, y_predict, average='macro')
            print("all_acc:", all_acc,"all_f1:", all_f1)

            labelsss = []
            y_predicts = []
            for i in labelss:
                labelsss.append(cell_types[i])
            for i in y_predict:
                y_predicts.append(cell_types[i])
            
            
            #Storing predicted cell types and the scTransformer
            log_dir = "log/"
            if (not os.path.isdir(log_dir)):
                os.makedirs(log_dir)
            
            np.save(log_dir  + 'y_tests.npy', labelsss)
            np.save(log_dir  + 'y_predicts.npy', y_predicts)
            torch.save(model.state_dict(), log_dir + 'scTransformer.tar')

            with open(log_dir + "resilt.txt", "a") as f:
                f.writelines("acc:" + str(all_acc) + "\n")
                f.writelines('f1:' + str(all_f1) + "\n")