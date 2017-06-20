import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import entropy
import sklearn.neighbors as NN

# read data

decayfile = "/net/mraid14/export/data/users/lubling/datasets/scell/results/esh/hyb_mm9_2i_is/tables/glob_decay_res_step0.125.txt"
replifile = "/net/mraid14/export/data/users/lubling/datasets/scell/results/esh/hyb_mm9_2i_is/tables/diploids_features_table.txt"

d = pd.read_csv(decayfile, sep="\t")
binn = d.shape[1]
d = d.rename(lambda e: e.split("scell.nextera.")[1])
r = pd.read_csv(replifile, sep="\t")
r = r.set_index('cell_nm')
df = d.join(r, how="inner")
assert np.all(df.columns[:d.shape[1]]==d.columns)
color = np.array(df['repli_score'])
no_score = np.isnan(df['repli_score'])
df = df[~no_score]
color = df['repli_score']

# extract features

y = np.array(df.iloc[:,range(1,binn-1)]).astype(np.float) # omit non-digest and trans
y[y==0] = 1 # add a pseudocount
y = y/y.sum(axis=1)[:,np.newaxis]

celln = y.shape[0]
binn = y.shape[1]

# compute distances

dists = np.zeros((celln, celln))
for i in range(celln):
    for j in range(i+1, celln):
        dists[i,j] = entropy(y[i,:], qk=y[j,:], base=None)
        dists[j,i] = entropy(y[j,:], qk=y[i,:], base=None)

symdist = dists+dists.T

# construct NN graph

k_nn = 7
nn = NN.NearestNeighbors(n_neighbors=k_nn, metric='precomputed', n_jobs=-1)
nn.fit(symdist)
dist2neighs, neighs =  nn.kneighbors()

adj = np.zeros((symdist.shape[0], symdist.shape[0]))
for i in range(symdist.shape[0]):
    for ki in range(k_nn):
        adj[i,neighs[i,ki]] = 1
adj = np.maximum(adj, adj.T)

# compute spectral embedding

lap = np.diag(adj.sum(axis=1))-adj
vals, vecs = np.linalg.eig(lap)
tmp = np.argsort(vals)
vals = vals[tmp]
vecs = vecs[:,tmp]

# plot

# fig = plt.figure(figsize=[2, 2])
plt.rcParams['xtick.labelsize'] = 8 
plt.rcParams['ytick.labelsize'] = 8 

fig = plt.figure()
tmporder = np.argsort(color)
plt.scatter(vecs[tmporder,1], vecs[tmporder,2], c=color[tmporder], marker='x', alpha=1, cmap=plt.cm.RdYlBu_r)
# plt.scatter(vecs[tmporder,1], vecs[tmporder,2], c=color[tmporder], marker='x', alpha=1, cmap=plt.cm.RdYlBu_r, s=8)
plt.colorbar()
plt.xlabel("axis 1", fontsize=10)
plt.ylabel("axis 2", fontsize=10)
plt.grid('on')
plt.xticks(plt.gca().get_xticks()[1::2])
fig.set_size_inches([5, 5],forward=True)
plt.savefig('fig4yaniv.png')

