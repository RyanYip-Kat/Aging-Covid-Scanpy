import scanpy as sc
import argparse
import sys
import os
sys.path.append("../")

from core.utils import get_rank_group_genes

parser=argparse.ArgumentParser()
parser.add_argument("-o","--outdir",type=str,default="output")
parser.add_argument("-d","--data",type=str,default=None,help="data")

args=parser.parse_args()

if not os.path.exists(args.outdir):
   os.makedirs(args.outdir)

adata=sc.read_h5ad(args.data)
sc.tl.rank_genes_groups(adata,groupby="status", method='t-test',n_genes=100,corr_method="bonferroni")
get_rank_group_genes(adata,args.outdir,50)
get_rank_group_genes(adata,args.outdir)
