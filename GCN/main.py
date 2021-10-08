import os
import time
import numpy as np
import pandas as pd 
import datatable as dt
from sklearn.metrics import accuracy_score, f1_score, roc_auc_score
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
import pytorch_lightning as pl
from pytorch_lightning.callbacks import LearningRateMonitor, ModelCheckpoint
os.chdir("/mnt/home/seguraab/Shiu_Lab/Collabs/Multi_Omic/multi-omics/GCN")
import simple_GCN

def main():

    geno_sub, X_train, train_loader, test_loader, test_tensor = simple_GCN.split_data()

    # Define the network
    net = simple_GCN.GCN(in_channels = geno_sub.shape[1], out_channels = [geno_sub.shape[1], 500, X_train.size()[0]])
    net = net.to(device="cuda:0")
    print(net)

    optimizer = torch.optim.SGD(net.parameters(), lr=0.2) # optimizer (variant of gradient descent)
    loss_func = torch.nn.MSELoss()  # mean squared loss function for regression

    start = time.time()
    # Train the network
    print("Training... ")
    for epoch in range(500):
        current_loss = 0.0
        print(f"\nTest: Epoch {epoch}")
        for input, batch in enumerate(train_loader, 0):
            optimizer.zero_grad()           # clear gradients for next train
            inputs, target = batch
            prediction = net(inputs)             # prediction based on input x
            loss = loss_func(prediction, target) # must be (1. nn output, 2. target)
            loss.backward()                 # backpropagation, compute gradients
            optimizer.step()                # apply gradients

            # print statistics
            current_loss += loss.item()
            print('Loss after mini-batch %5d: %.3f' % (input + 1, current_loss / 500))
            current_loss = 0.0
    
    end = time.time()   
    elapsed = end - start
    print("Training complete; Elapsed time: ", elapsed)

    # Evaluate trained model on test set
    

if __name__ == "__main__":
    main()