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
import pytorch_lightning as pl
from pytorch_lightning.callbacks import LearningRateMonitor, ModelCheckpoint
os.chdir("/mnt/home/seguraab/Shiu_Lab/Collabs/Multi_Omic/multi-omics/GCN")
import simple_GCN

geno_data="SNP_binary_matrix_383_accessions_drop_all_zero_MAF_larger_than_0.05_converted.csv"
pheno_data="Phenotype_value_383_common_accessions_2017_Grimm.csv"

geno = dt.fread(geno_data) # read in genotype data
geno = geno.to_pandas() # convert dataframe to pandas dataframe
geno = geno.sort_values(by=geno.columns[0], axis=0) # sort values by sample ID
geno = geno.set_index(geno.columns[0], drop=True) # set index to sample ID
geno_sub = geno.iloc[:,0:100000]
features = geno_sub.columns # columns as features

pheno = pd.read_csv(pheno_data, index_col=0) # read in phenotype data
label = pheno.FT10_mean # flowering time as label

# Split geno and pheno into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(geno_sub, label, test_size=64)

# Convert to PyTorch tensors
X_train = torch.tensor(X_train.values.astype(np.float32))
X_test = torch.tensor(X_test.values.astype(np.float32))
y_train = torch.tensor(y_train.values.astype(np.float32))
y_test = torch.tensor(y_test.values.astype(np.float32))
#X_train.size()
train_tensor = data.TensorDataset(X_train, y_train) # Tensor training dataset to feed
train_loader = data.DataLoader(dataset = train_tensor, batch_size = 10, shuffle = True) # data loader to feed training data to network

# Thought: Do I need to normalize my data?

def main():
    
    # Define the network
    net = simple_GCN.GCN(in_channels = geno.shape[1], out_channels = [geno_sub.shape[1], 50000, X_train.size()[0]])
    print(net)

    optimizer = torch.optim.SGD(net.parameters(), lr=0.2) # optimizer (variant of gradient descent)
    loss_func = torch.nn.MSELoss()  # mean squared loss function for regression

    # Train the network
    print("Training... ")
    for epoch in range(500):
        print("\nTest: Epoch {:d}".format(epoch))
        for batch in train_loader:
            optimizer.zero_grad()           # clear gradients for next train
            input, target = batch
            prediction = net(x)             # prediction based on input x
            loss = loss_func(prediction, y) # must be (1. nn output, 2. target)
            loss.backward()                 # backpropagation, compute gradients
            optimizer.step()                # apply gradients
        

    # Apply trained model to test set

if __name__ == "__main__":
    main()