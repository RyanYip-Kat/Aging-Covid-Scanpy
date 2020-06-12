import argparse
import os
import scanpy as sc
import numpy as np
import pandas as pd
import sys

sys.path.append("../")

from core.utils import subset_by_column,subset_by_cell_feature
from collections import Counter

parser=argparse.ArgumentParser()
parser.add_argument("-o","--outdir",type=str,default=None)
parser.add_argument("-d","--data",type=str,default=None)
parser.add_argument("-s","--subset",nargs="+",type=str,default=None)
parser.add_argument("-c","--column",type=str,default="leiden")
parser.add_argument("-i","--invert", action='store_true', default=False)

args=parser.parse_args()
if not os.path.exists(args.outdir):
   os.makedirs(args.outdir)

sc.settings.figdir=args.outdir
print("*** Loading data")
adata=sc.read_h5ad(args.data)

print("Subset data by {}".format(args.column))
if args.column is not None and args.subset is not None:
    print("*** subset")
    print("Subset data by {}".format(args.column))
    adata=subset_by_column(adata,args.subset,args.column,invert=args.invert)
    
X=adata.raw.X.toarray()
matrix=pd.DataFrame(X.transpose(),columns=adata.obs_names,index=adata.var_names)

filename=os.path.join(args.outdir,"expression_matrix.txt")
print("*** Expression matrix output in {}".format(filename))
matrix.to_csv(filename,sep="\t")

print("*** Done! ***")
