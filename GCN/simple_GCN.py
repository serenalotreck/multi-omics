import os
import time
import numpy as np
import pandas as pd 
import datatable as dt
from sklearn.metrics import accuracy_score, f1_score, roc_auc_score
from sklearn.model_selection import train_test_split
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.utils.data as data
import torch.optim as optim
from torch_geometric.nn import GCNConv
import pytorch_lightning as pl
from pytorch_lightning.callbacks import LearningRateMonitor, ModelCheckpoint

# Set seed for reproducibility
#pl.seed_everything(42)

# Ensure that all operations are deterministic on GPU for reproducibility
torch.backends.cudnn.determinstic = True
torch.backends.cudnn.benchmark = False

device = torch.device("cuda:0") if torch.cuda.is_available() else torch.device("cpu") # set to gpu or cpu
print(device)
torch.cuda.device_count() # number of gpus

# Set path to directory containing data
PATH = "/mnt/home/seguraab/Shiu_Lab/Collabs/Multi_Omic"
os.chdir(PATH)

# Define graph convolutional network
class GCN(nn.Module):
    def __init__(self, in_channels, out_channels): # size of input, size of outputs
        super(GCN, self).__init__()
        torch.manual_seed(42) # for reproducibility
        # Add layers
        self.layer1 = GCNConv(in_channels, out_channels[0], cached = True)
        self.layer2 = GCNConv(out_channels[0], out_channels[1], cached = True)
        self.layer3 = GCNConv(out_channels[1], out_channels[2], cached = True)
        self.classifier = nn.Sequential(nn.Linear(in_channels, out_channels[2]))

    def forward(self, x):
        x = self.layer1(x) # graph convolutional layer
        x = F.leaky_relu(x, 0.25) # activation function
        x = F.dropout(x, self.dropout, training = self.training) # regularization
        x = self.layer2(x)
        x = F.leaky_relu(x, 0.25)
        x = F.dropout(x, self.dropout, training=self.training)
        x = self.layer3(x)
        x = F.leaky_relu(x, 0.25)

        # Apply the linear classifier
        out = self.classifier(x)

        return out, x   
 
 # include function to split data

def split_data(): # geno_data, pheno_data args; move to simple_GCN
    geno_data="SNP_binary_matrix_383_accessions_drop_all_zero_MAF_larger_than_0.05_converted.csv"
    pheno_data="Phenotype_value_383_common_accessions_2017_Grimm.csv"

    geno = dt.fread(geno_data) # read in genotype data
    geno = geno.to_pandas() # convert dataframe to pandas dataframe
    geno = geno.sort_values(by=geno.columns[0], axis=0) # sort values by sample ID
    geno = geno.set_index(geno.columns[0], drop=True) # set index to sample ID
    geno_sub = geno.iloc[:,0:1000]
    features = geno_sub.columns # columns as features

    pheno = pd.read_csv(pheno_data, index_col=0) # read in phenotype data
    label = pheno.FT10_mean # flowering time as label

    # Split geno and pheno into training and testing sets
    X_train, X_test, y_train, y_test = train_test_split(geno_sub, label, test_size=64)

    # Convert to PyTorch tensors
    X_train = torch.tensor(X_train.values.astype(np.float32), device = "cuda:0")
    X_test = torch.tensor(X_test.values.astype(np.float32), device = "cuda:1")
    y_train = torch.tensor(y_train.values.astype(np.float32), device = "cuda:0")
    y_test = torch.tensor(y_test.values.astype(np.float32), device = "cuda:1")
    
    train_tensor = data.TensorDataset(X_train, y_train) # Tensor training dataset to feed
    train_loader = data.DataLoader(dataset = train_tensor, batch_size = 10, shuffle = True, pin_memory = False) # data loader to feed training data to network
    test_tensor = data.TensorDataset(X_test, y_test)
    test_loader = data.DataLoader(dataset=test_tensor, shuffle = True, pin_memory = False)
    return geno_sub, X_train, train_loader, test_loader, test_tensor # need to figure what to return

# Thought: Do I need to normalize my data?
