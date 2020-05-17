import scanpy as sc
import argparse
import sys
import os
sys.path.append("../")

from core.utils import get_rank_group_genes,subset_by_column

parser=argparse.ArgumentParser()
parser.add_argument("-o","--outdir",type=str,default="output")
parser.add_argument("-d","--data",type=str,default=None,help="data")
parser.add_argument("-s","--subset",nargs="+",type=str,default=None,help="subset")
parser.add_argument("-c","--column",type=str,default="leiden",help="data")

args=parser.parse_args()

if not os.path.exists(args.outdir):
   os.makedirs(args.outdir)

adata=sc.read_h5ad(args.data)
if args.subset is not None:
    adata=subset_by_column(adata,args.subset,args.column)
sc.tl.rank_genes_groups(adata,groupby="status", method='t-test',n_genes=500,corr_method="bonferroni")
get_rank_group_genes(adata,args.outdir,50)
get_rank_group_genes(adata,args.outdir)
