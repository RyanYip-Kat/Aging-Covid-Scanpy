import scanpy as sc
import os
import pandas as pd
import argparse
from collections import Counter

parser=argparse.ArgumentParser()
parser.add_argument("-o","--outdir",type=str,default="output")
parser.add_argument("-c","--column",type=str,default="leiden")
args=parser.parse_args()

if not os.path.exists(args.outdir):
   os.makedirs(args.outdir)

adata=sc.read_h5ad("../20200509/output/aging-11/BC/adata.h5ad")
DATA=pd.read_csv("MEMORY.csv",sep=",")
Barcodes=DATA.Barcode.values

metadata=adata.obs
metadata["barcode"]=metadata.index
clusters=metadata[["barcode",args.column]]

print("------------")
print(Counter(metadata[args.column]))

new_Clusters=[]
for barcode,column in clusters.values:
    if barcode in Barcodes:
        new_Clusters.append("Memory")
    else:
        new_Clusters.append(column)
print("-------------")
print(Counter(new_Clusters))

adata.obs["new_Clusters"]=new_Clusters
adata.write(os.path.join(args.outdir,'adata.h5ad'), compression='gzip')
