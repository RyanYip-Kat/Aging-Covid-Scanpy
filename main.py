import scanpy as sc
from core.core import Model
from core.ScanpyConfig import *

path="/home/ye/Work/BioAligment/10X/output/20200407_aging-XXX"
net=Model(path)
print(net)

adata=net._load_data()
adata=net._subset(adata,[10,13])
adata=net._process(adata)
adata=net._reduction(adata)
adata=net._FindMarkers(adata)

adata.write('./data/adata.h5ad', compression='gzip')
