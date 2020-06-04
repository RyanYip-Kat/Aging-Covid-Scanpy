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
parser.add_argument("-g","--groupby",type=str,default="status",help="data")

args=parser.parse_args()

if not os.path.exists(args.outdir):
   os.makedirs(args.outdir)

adata=sc.read_h5ad(args.data)
if args.subset is not None:
    adata=subset_by_column(adata,args.subset,args.column)
    adata=adata.raw.to_adata()
    adata.raw=adata
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.highly_variable_genes(adata,n_top_genes=5000)
    adata= adata[:, adata.var.highly_variable]
    sc.pp.scale(adata, max_value=10)

sc.tl.rank_genes_groups(adata,groupby=args.groupby,method="wilcoxon",n_genes=100,corr_method="bonferroni")
get_rank_group_genes(adata,pval=None,fc=None,outdir=args.outdir,n_top=50)
get_rank_group_genes(adata,pval=None,fc=None,outdir=args.outdir)