import time
import numpy as np
import pandas as pd 
import datatable as dt
from sklearn.metrics import accuracy_score, f1_score, roc_auc_score
from sklearn.model_selection import train_test_split
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
import pytorch_lightning as pl
from pytorch_lightning.callbacks import LearningRateMonitor, ModelCheckpoint

# Set seed for reproducibility
pl.seed_everything(42)

# Ensure that all operations are deterministic on GPU for reproducibility
torch.backends.cudnn.determinstic = True
torch.backends.cudnn.benchmark = False

device = torch.device("cuda:0") if torch.cuda.is_available() else torch.device("cpu")
print(device)

# Path to directory containing data
PATH = "/mnt/home/seguraab/Shiu_Lab/Collabs/Multi_Omic"

# Define graph convolutional network
class GCN(nn.Module):
    def __init__(self, in_channels, out_channels): # size of input, size of outputs
        super(GCN, self).__init__()
        torch.manual_seed(42) # for reproducibility
        # Add layers
        self.layer1 = GCNConv(in_channels, out_channels[0])
        self.layer2 = GCNConv(out_channels[0], out_channels[1])
        self.layer3 = GCNConv(out_channels[1], out_channels[2])
        self.classifier = nn.Sequential(Linear(in_channels, out_channels[2]))

    def forward(self, x, adj_mat):
        x = self.layer1(x, adj_mat) # graph convolutional layer
        x = F.leaky_relu(x, 0.25) # activation function
        x = F.dropout(x, self.dropout, training = self.training) # regularization
        x = self.layer2(x, adj_mat)
        x = F.leaky_relu(x, 0.25)
        x = F.dropout(x, self.dropout, training=self.training)
        x = self.layer3(x, adj_mat)
        x = F.leaky_relu(x, 0.25)

        # Apply the linear classifier
        out = self.classifier(x)

        return out, x   


def split_train_test(geno_data, pheno_data):
    geno = dt.fread(geno_data) # read in genotype data
    geno = geno.to_pandas() # convert dataframe to pandas dataframe
    features = geno.columns[1:] # columns as features

    pheno = pd.read_csv(pheno_data) # read in phenotype data
    geno = geno.rename(columns={geno.columns[0]:pheno.columns[0]}) # ensure first columns match

    label = pheno.FT10_mean # flowering time as label

    # Split geno and pheno into training and testing sets
    X_train, X_test, y_train, y_test = train_test_split(geno, label, test_size=64)

    # Convert to PyTorch tensors
    X_train = torch.FloatTensor(X_train)
    X_test = torch.FloatTensor(X_test)
    y_train = torch.FloatTensor(y_train)
    y_test = torch.FloatTensor(y_train)

    return X_train, X_test, y_train, y_test

def main():
    
    # Define the network
    net = GCN(in_channels = , out_channels = )
    print(net)

    optimizer = torch.optim.SGD(net.parameters(), lr=0.2) # optimizer (variant of gradient descent)
    loss_func = torch.nn.MSELoss()  # mean squared loss function for regression

    # Train the network
    print("Training... ")
    prediction = net(x) # prediction based on input x
    
    loss = loss_func(prediction, y) # must be (1. nn output, 2. target)
    
    optimizer.zero_grad()   # clear gradients for next train
    loss.backward()         # backpropagation, compute gradients
    optimizer.step()        # apply gradients

    # Apply the model to the Test set
    for epoch in range(500):
        print("\nTest: Epoch {:d}".format(epoch))
        
        
        print("Test ACC: {:.3f}".format(accuracy_score(labels_trte[trte_idx["te"]], te_prob.argmax(1))))
        print("Test F1: {:.3f}".format(f1_score(labels_trte[trte_idx["te"]], te_prob.argmax(1))))
        print("Test AUC: {:.3f}".format(roc_auc_score(labels_trte[trte_idx["te"]], te_prob[:,1])))
        



if __name__ == "__main__":
    main()

