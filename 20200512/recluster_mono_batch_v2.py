import argparse
import os
import scanpy as sc
import numpy as np
import sys

sys.path.append("../")

from core.tl import reduction,process
from core.core import Model
from core.utils import subset_by_column,subset_by_cell_feature

parser=argparse.ArgumentParser()
parser.add_argument("-o","--outdir",type=str,default=None)
parser.add_argument("-d","--data",type=str,default=None)
parser.add_argument("-s","--subset",nargs="+",type=str,default=None)
parser.add_argument("-c","--column",type=str,default="leiden")

args=parser.parse_args()
if not os.path.exists(args.outdir):
   os.makedirs(args.outdir)

print("*** Loading data")
adata=sc.read_h5ad(args.data)

print("*** subset")
print("Subset data by {}".format(args.column))
adata=subset_by_column(adata,args.subset,args.column)
adata=adata.raw.to_adata()
idents=adata.obs["idents"]
batch=[]
for ident in idents.values:
    if ident=="2":
        batch.append("CAT-4")
    elif ident=="3":
        batch.append("PCV-1")
    elif ident=="4":
        batch.append("PCV-2")
    elif ident=="7":
        batch.append("NPDR-1")
    elif ident=="8":
        batch.append("WAMD-1")
    elif ident=="5":
        batch.append("PCV-3")
    else:
        batch.append("Others")

adata.obs["batch"]=batch

print("*** tranform and reduction")
adata=process(adata,key="batch")
adata=reduction(adata)

print("*** find markers")
sc.tl.rank_genes_groups(adata,groupby="leiden", method='t-test',n_genes=500,corr_method="bonferroni")
adata.write(os.path.join(args.outdir,'adata.h5ad'), compression='gzip')
