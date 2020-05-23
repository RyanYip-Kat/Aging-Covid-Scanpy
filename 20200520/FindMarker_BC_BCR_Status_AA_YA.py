import scanpy as sc
import argparse
import sys
import os
import pandas as pd
import numpy as np
sys.path.append("../")
from core.utils import get_rank_group_genes,subset_by_column
from collections import Counter
parser=argparse.ArgumentParser()
parser.add_argument("-o","--outdir",type=str,default="output")

args=parser.parse_args()

if not os.path.exists(args.outdir):
   os.makedirs(args.outdir)

sc.settings.figdir=args.outdir
adata=sc.read_h5ad("../20200513/output/aging-11/Total/adata_doublet.h5ad")
DATA=pd.read_csv("BCR-DF.csv",sep=",")

status=["AA","YA"]
for s in status:
    adata_subset=adata[adata.obs["status"]==s]
    df=DATA[DATA["status"]==s]
    barcode=df["Barcode"].values.tolist()
    celltype=[]
    for b in adata_subset.obs_names:
        if b in barcode:
            celltype.append("BCR")
        else:
            celltype.append("NOT_BCR")

    print(Counter(celltype))
    adata_subset.obs["BCR"]=celltype
    sc.tl.rank_genes_groups(adata_subset,groupby="BCR", method='t-test',n_genes=500,corr_method="bonferroni")
    sc.tl.filter_rank_genes_groups(adata_subset,groupby="BCR",min_fold_change=0.25)

    get_rank_group_genes(adata_subset,args.outdir,prefix=s)
    
