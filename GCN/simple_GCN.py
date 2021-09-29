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
import torch.optim as optim
from torch_geometric.nn import GCNConv
import pytorch_lightning as pl
from pytorch_lightning.callbacks import LearningRateMonitor, ModelCheckpoint

# Set seed for reproducibility
#pl.seed_everything(42)

# Ensure that all operations are deterministic on GPU for reproducibility
torch.backends.cudnn.determinstic = True
torch.backends.cudnn.benchmark = False

device = torch.device("cuda:0") if torch.cuda.is_available() else torch.device("gpu") # set to gpu or cpu
print(device)

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
 
 # include function to split data
