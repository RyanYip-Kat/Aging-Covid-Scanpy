import scanpy as sc
import os
import pandas as pd
import argparse

parser=argparse.ArgumentParser()
#parser.add_argument("-o","--outdir",type=str,default="output")
parser.add_argument("-d","--data",type=str,default=None,help="data")
parser.add_argument("-f","--barcode",type=str,default=None,help="the csv file of barcode to be specified")
parser.add_argument("-c","--column",type=str,default="celltype",help="the columns in adata.obs to be used")
parser.add_argument("-n","--name",type=str,default=None)

args=parser.parse_args()

#if not os.path.exists(args.outdir):
#   os.makedirs(args.outdir)

print("*** Loading data from : {}".format(args.data))
adata=sc.read_h5ad(args.data)
DATA=pd.read_csv(args.barcode,sep=",")

barcodes=DATA["Barcode"].values.tolist()
metadata=adata.obs
metadata["cells"]=metadata.index

new_celltype=[]
for barcode,celltype in zip(metadata["cells"].values.tolist(),metadata[args.column].values.tolist()):
    if barcode in barcodes:
        new_celltype.append(args.name)
    else:
        new_celltype.append(celltype)

adata.obs["new_celltype"]=new_celltype
adata.write(args.data)


