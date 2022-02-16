""" 
GCN for Classification from praxidike97 on GitHub adapted for regression with genotype data.

GCN source: https://github.com/praxidike97/GraphNeuralNet/blob/master/main.py
layer source: https://github.com/dizhu-gis/SRGCNN/blob/main/SRGCNN_demo.ipynb
Adjacency matrix source: https://github.com/txWang/MOGONET

To run on command line if CUDA error (device-side assert triggered): 
    CUDA_LAUNCH_BLOCKING=1 python praxidike97_gcn_regression_changed.py 

Data used:
    Planetoid: 
        Data(x=[2708, 1433], edge_index=[2, 10556], y=[2708], train_mask=[2708], val_mask=[2708], test_mask=[2708])
        num_node_features = 1433; num_classes = 7
    
    Multi-omics genotype: Data(x=[383, 1000], edge_index=[2, 5744], y=[383], train_mask=[383], test_mask=[383])
        1000 features, 383 instances
    
    Yeast genotype: Data(x=[750, 64456], edge_index=[2, 1251], y=[750], train_mask=[750], test_mask=[750])
        64456 features, 750 instances
"""

from pyexpat.errors import XML_ERROR_INVALID_TOKEN
import torch
from torch import nn
from torch import tensor
import torch.nn.functional as F
from torch_geometric.datasets import Planetoid
from torch_geometric.nn import MessagePassing
from torch_geometric.utils import add_self_loops, degree
from torch_geometric.data import Data
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import time
import pandas as pd
import datatable as dt

load_from_preprocess = False

# Graph Convolutional layer
class GCNConv(MessagePassing):
    def __init__(self, in_channels, out_channels):
        super(GCNConv, self).__init__(aggr='add')  # "Add" aggregation
        self.lin = torch.nn.Linear(in_channels, out_channels)

    def forward(self, x, edge_index):
        # Step 1: Add self-loops
        edge_index, _ = add_self_loops(edge_index, num_nodes=x.size(0))

        # Step 2: Multiply with weights
        x = self.lin(x)

        # Step 3: Calculate the normalization
        row, col = edge_index
        deg = degree(row, x.size(0), dtype=x.dtype)
        deg_inv_sqrt = deg.pow(-0.5)
        norm = deg_inv_sqrt[row] * deg_inv_sqrt[col]

        # Step 4: Propagate the embeddings to the next layer
        return self.propagate(edge_index, size=(x.size(0), x.size(0)), x=x,
                              norm=norm)

    def message(self, x_j, norm):
        # Normalize node features.
        return norm.view(-1, 1) * x_j

# Graph convolutional layer for regression
class GCNConv_reg(nn.Module):
    """Basic graph convolution operation that incorporate both spatial lagged X 
    and spatial lagged Y (to be used in the basic SRGCNNs model)"""
    def __init__(self, f_in, f_out, use_bias=True, activation=F.relu):#hidden layer with relu activation
        super().__init__()
        self.f_in = f_in
        self.f_out = f_out
        self.use_bias = use_bias
        self.activation = activation
        self.weight = nn.Parameter(torch.FloatTensor(f_in, f_out),requires_grad=True)###requires_grad:whether allow weights to be updated
        self.bias = nn.Parameter(torch.FloatTensor(f_out)) if use_bias else None
        self.initialize_weights()
    
    def initialize_weights(self):
        nn.init.constant_(self.weight,1)
        if self.use_bias: nn.init.constant_(self.bias,0)
        
    def forward(self, input, adj):
        support = torch.mm(input, self.weight)
        output = torch.mm(adj, support) #adj here has to be renormalized     
        
        if self.use_bias: output.add_(self.bias)
        if self.activation is not None: output=self.activation(output) 
        
        return output

# GCN Model
class NetReg(torch.nn.Module):
    def __init__(self, num_node_features):
        super().__init__()
        # for classification
        self.conv1 = GCNConv(data.x.shape[1], 200) # data.x.shape[1] = 750 samples; dataset.num_node_features, 16
        self.conv2 = GCNConv(200, 1) # 16, 1
        # for regression
        #self.conv1 = GCNConv(1000, 400) 
        #self.conv2 = GCNConv(400, 383)

    def forward(self, data):
        x, edge_index = data.x, data.edge_index
        x = self.conv1(x, edge_index)
        x = F.relu(x)
        x = F.dropout(x, training=self.training)
        #print('dropout x', x, x.size())
        x = self.conv2(x, edge_index)
        #print('x 2', x, x.size())
        
        #return x # for regression
        return x# for classification
        

class Net(torch.nn.Module): # Classification
    def __init__(self, dataset):
        super(Net, self).__init__()
        # for classification
        self.conv1 = GCNConv(dataset.num_node_features, 16) 
        self.conv2 = GCNConv(16, dataset.num_classes)


    def forward(self, data):
        x, edge_index = data.x, data.edge_index
        x = self.conv1(x, edge_index)
        x = F.relu(x) # this is messing up the structure of the tensor for genotype data!! :( but it's fine for planetoid. Why is that?
        x = F.dropout(x, training=self.training)
        x = self.conv2(x, edge_index)
        return F.log_softmax(x, dim=1) # for classification

# Graph of dataset
def plot_dataset(dataset):
    edges_raw = dataset.data.edge_index.numpy()
    edges = [(x, y) for x, y in zip(edges_raw[0, :], edges_raw[1, :])]
    labels = dataset.data.y.numpy()

    G = nx.Graph()
    G.add_nodes_from(list(range(np.max(edges_raw))))
    G.add_edges_from(edges)
    plt.subplot(111)
    options = {
                'node_size': 30,
                'width': 0.2,
    }
    nx.draw(G, with_labels=False, node_color=labels.tolist(), cmap=plt.cm.tab10, font_weight='bold', **options)
    plt.savefig("planetoid_graph.png")


# Evaluate the model on the test set
def test(data, train=True, val=False):
    model.eval()
    correct = 0
    pred = model(data)#.max(dim=1)[1]
    if train:
        correct += pred[data.train_mask].eq(data.y[data.train_mask]).sum().item()
        return correct / (len(data.y[data.train_mask]))
    if val:
        correct += pred[data.val_mask].eq(data.y[data.val_mask]).sum().item()
        return correct / (len(data.y[data.val_mask]))
    else:
        correct += pred[data.test_mask].eq(data.y[data.test_mask]).sum().item()
        return correct / (len(data.y[data.test_mask]))

# Evaluate the model on the validation or test set
def test_reg(data, train=True):
    model.eval()
    correct = 0
    pred = model(data)
    return (pred - data.y).pow(2).sum().sqrt().item()

# train the model
def train_reg(data, plot=False):    
    train_accuracies, val_accuracies = list(), list()
    start = time.time()
    for epoch in range(100):
        model.train()
        optimizer.zero_grad()
        out = model(data)
        
        # for regression
        #loss = torch.nn.MSELoss()
        #output = loss(out[data.train_mask], data.y[data.train_mask]) 
        #output.backward()

        # for regression
        loss = F.mse_loss(out[data.val_mask],data.y[data.val_mask])
        loss.backward()
        optimizer.step()

        train_acc = test_reg(data)
        val_acc = test_reg(data, train=False)

        train_accuracies.append(train_acc)
        val_accuracies.append(val_acc)
        print('Epoch: {:03d}, Loss: {:.5f}, Train Acc: {:.5f}, Val Acc: {:.5f}'.
              format(epoch, loss, train_acc, val_acc))
    end = time.time()
    print("Elapsed time: ", end-start)

    # Test accuracy
    test_acc = test(data, train=False)
    print("Test Acc: {:.5f}".format(test_acc))

    if plot: # plot AUC curve
        plt.plot(train_accuracies, label="Train mse loss")
        plt.plot(val_accuracies, label="Validation mse loss")
        plt.ylim([0,1000])
        plt.xlabel("# Epoch")
        plt.ylabel("mse")
        plt.title("new Dataset")
        plt.legend(loc='upper right')
        plt.savefig("auc_praxidike_gcn_yeast_val.png")

# train the model
def train(data, plot=False):    
    train_accuracies, val_accuracies = list(), list()
    start = time.time()
    for epoch in range(100):
        model.train()
        optimizer.zero_grad()
        out = model(data)

        # for classification
        loss = F.nll_loss(out[data.val_mask], data.y[data.val_mask])
        loss.backward()
        optimizer.step()

        train_acc = test(data)
        val_acc = test(data, train=False, val=True)

        train_accuracies.append(train_acc)
        val_accuracies.append(val_acc)
        print('Epoch: {:03d}, Loss: {:.5f}, Train Acc: {:.5f}, Val Acc: {:.5f}'.
              format(epoch, loss, train_acc, val_acc))
   
    end = time.time()
    print("Elapsed time: ", end-start)

    # Test accuracy
    test_acc = test(data, train=False)
    print("Test Acc: {:.5f}".format(test_acc))

    if plot: # plot AUC curve
        plt.plot(train_accuracies, label="Train accuracy")
        plt.plot(val_accuracies, label="Validation accuracy")
        plt.xlabel("# Epoch")
        plt.ylabel("Accuracy")
        plt.title("Planetoid Dataset")
        plt.legend(loc='upper right')
        plt.savefig("auc_praxidike_gcn_val.png")

'''
MOGONET Method for Computing Adjacency Matrix
1. Calculate adjacency matrix parameter (cal_adj_mat_parameter)
    a. Compute cosine similarity matrix (cosine_distance_torch)
2. Generate adjacency matrix tensor (gen_adj_mat_tensor)
    a. Compute cosine similarity matrix (cosine_distance_torch)
    b. Generate graph from distance matrix (graph_from_dist_torch)
    c. Compute identity matrix
    d. Compute adjacency matrix (not equation 2 in paper)
3. Generate test adjacency matrix tensor (gen_test_adj_mat_tensor)
'''
def cosine_distance_torch(x1, x2=None, eps=1e-8):
    x2 = x1 if x2 is None else x2
    w1 = x1.norm(p=2, dim=1, keepdim=True)
    w2 = w1 if x2 is x1 else x2.norm(p=2, dim=1, keepdim=True)
    return 1 - torch.mm(x1, x2.t()) / (w1 * w2.t()).clamp(min=eps)

def cal_adj_mat_parameter(edge_per_node, data, metric="cosine"):
    assert metric == "cosine", "Only cosine distance implemented"
    dist = cosine_distance_torch(data, data)
    parameter = torch.sort(dist.reshape(-1,)).values[edge_per_node*data.shape[0]]
    return parameter.data.cpu().numpy().item()

def graph_from_dist_tensor(dist, parameter, self_dist=True):
    if self_dist:
        assert dist.shape[0]==dist.shape[1], "Input is not pairwise dist matrix"
    g = (dist <= parameter).float() # binary 0s and 1s for False/True
    if self_dist:
        diag_idx = np.diag_indices(g.shape[0])
        g[diag_idx[0], diag_idx[1]] = 0
    return g

def to_sparse(x):
    x_typename = torch.typename(x).split('.')[-1]
    sparse_tensortype = getattr(torch.sparse, x_typename)
    indices = torch.nonzero(x)
    if len(indices.shape) == 0:  # if all elements are zeros
        return sparse_tensortype(*x.shape)
    indices = indices.t()
    values = x[tuple(indices[i] for i in range(indices.shape[0]))]
    return sparse_tensortype(indices, values, x.size()), indices, values

def gen_adj_mat_tensor(data, parameter, metric="cosine"):
    assert metric == "cosine", "Only cosine distance implemented"
    dist = cosine_distance_torch(data, data)
    g = graph_from_dist_tensor(dist, parameter, self_dist=True)
    if metric == "cosine":
        adj = 1-dist
    else:
        raise NotImplementedError
    adj = adj*g 
    adj_T = adj.transpose(0,1)
    I = torch.eye(adj.shape[0])
    cuda = True if torch.cuda.is_available() else False
    adj = adj + adj_T*(adj_T > adj).float() - adj*(adj_T > adj).float()
    adj = F.normalize(adj + I, p=1)
    adj = to_sparse(adj)   
    return adj



# create PyTorch Data object using genotype data
def test_geno():
    # Datasets
    #path = "/mnt/home/seguraab/Shiu_Lab/Collabs/Multi_Omic/Data" # multi-omics dataset
    #geno_data="%s/SNP_binary_matrix_383_accessions_drop_all_zero_MAF_larger_than_0.05_converted.csv"%path
    #pheno_data="%s/Phenotype_value_383_common_accessions_2017_Grimm.csv"%path
    #test_mask="%s/test_20perc.txt"%path
    path = "/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018"
    geno_data = "%s/geno.csv"%path
    pheno_data = "%s/pheno.csv"%path
    test_mask="%s/Test.txt"%path
    test_mask = pd.read_csv(test_mask, header=None)
    
    if not load_from_preprocess:
        geno = dt.fread(geno_data) # read in genotype data
        geno = geno.to_pandas() # convert dataframe to pandas dataframe  # 383, 1771291 
        # geno = pd.read_csv(geno_data)
        geno = geno.sort_values(by=geno.columns[0], axis=0) # sort values by sample ID
        geno = geno.set_index(geno.columns[0], drop=True) # set index to sample ID
        geno_sub = geno
        features = geno_sub.columns # columns as features
        
        pheno = pd.read_csv(pheno_data, index_col=0) # read in phenotype data
        #label = pheno.FT10_mean
        label = pheno.YPACETATE
        
        # Split geno and pheno into training and testing sets
        X_train = geno_sub.loc[~geno_sub.index.isin(test_mask[0])]
        X_test = geno_sub.loc[geno_sub.index.intersection(test_mask[0])]
        y_train = pheno.loc[~pheno.index.isin(test_mask[0])]
        y_test = pheno.loc[pheno.index.intersection(test_mask[0])]

        # Split training sets into training and validation sets
        X_val = X_train.sample(frac=0.2, random_state=25)
        X_train = X_train.loc[~X_train.index.isin(X_val.index)]
        y_val = y_train.loc[y_train.index.isin(X_val.index)]
        y_train = y_train.loc[~y_train.index.isin(y_val.index)]
        
        # Create masks
        train_mask = tensor([i in np.array(X_train.index) for i in np.array(geno.index)])
        val_mask = tensor([i in np.array(X_val.index) for i in np.array(geno.index)])
        test_mask = tensor([i in np.array(X_test.index) for i in np.array(geno.index)])
        
        # Convert to PyTorch tensors
        geno_sub = torch.tensor(geno_sub.values.astype(np.float32))
        #geno_sub[geno_sub==-1] = 0 # convert -1 to 0, just in case
        label = torch.tensor(label.values.astype(np.float32))
        X_train = torch.tensor(X_train.values.astype(np.float32))
        X_val = torch.tensor(X_val.values.astype(np.float32))
        X_test = torch.tensor(X_test.values.astype(np.float32))
        y_train = torch.tensor(y_train.YPACETATE.values) #FT10_mean
        y_val = torch.tensor(y_val.values.astype(np.float32))
        y_test = torch.tensor(y_test.YPACETATE.values) #FT10_mean
        
    else:
        geno_sub = torch.load(open("geno_sub.pth",'rb'))
        label = torch.load(open("label.pkl",'rb'))
        train_mask = torch.load(open("train_mask.pkl",'rb'))
        test_mask  = torch.load(open("test_mask.pkl",'rb'))
        X_train = torch.load(open("x_train.pkl",'rb'))
        X_test  = torch.load(open("x_test.pkl",'rb'))
        y_test  = torch.load(open("y_test.pkl",'rb'))
        y_train = torch.load(open("y_train.pkl",'rb'))

    # Compute adjacency matrix
    adj_parameter = 2 # edge_per_node
    adj_parameter_adaptive = cal_adj_mat_parameter(adj_parameter, X_train, "cosine")
    adj_train, edge_index, adj_values = gen_adj_mat_tensor(X_train, adj_parameter_adaptive, "cosine")

    # Create PyTorch geometric Data object
    return Data(x=geno_sub, edge_index=edge_index, y=label, train_mask=train_mask, val_mask=val_mask, test_mask=test_mask)

if __name__ == "__main__":
    import gc
    gc.collect()
    torch.cuda.empty_cache()
    device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu') # check device
    print(device)

    # planetoid data for classification
    # dataset = Planetoid(root='/tmp/Cora', name='Cora') # load dataset (2708 nodes/input tensor vars, 1433 instances)
    # data = dataset[0].to(device)
    
    # genotype data for regression
    dataset = test_geno()
    data = dataset.cuda(device)
    print("data", data)
    # model = Net(dataset) # create GCN model
    # model.cuda(device) # send to gpu
    model = NetReg(data.x.shape[1])
    model.cuda(device) # send to gpu
    # # Optimizer    
    optimizer = torch.optim.Adam(model.parameters(), lr=0.01, weight_decay=5e-4)
    # # Train the model
    train(data,plot=True)
    #train_reg(data, plot=True)