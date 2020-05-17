import scanpy as sc
import argparse
import sys
import os
sys.path.append("../")

from core.utils import get_rank_group_genes,subset_by_column

parser=argparse.ArgumentParser()
parser.add_argument("-o","--outdir",type=str,default="output")
parser.add_argument("-d","--data",type=str,default=None,help="data")

args=parser.parse_args()

if not os.path.exists(args.outdir):
   os.makedirs(args.outdir)

adata=sc.read_h5ad(args.data)

subset_lists={"P1":["0"],
              "P2":["1"],
              "P3":["3"],
              "P4":["4"],
              "P5":["0","1"],
              "P6":["3","4","6"],
              "P7":["6"],
              "P8":["10"],
              "P9":["3","4","6","10"]
              }

col_use="leiden"
for key,subset in subset_lists.items(): 
    adata_subset=subset_by_column(adata,subset,col_use,False)
    sc.tl.rank_genes_groups(adata_subset,groupby="status", method='t-test_overestim_var',n_genes=500,corr_method="bonferroni")
    
    outdir=os.path.join(args.outdir,key)
    
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    get_rank_group_genes(adata_subset,outdir,50)
    get_rank_group_genes(adata_subset,outdir)
