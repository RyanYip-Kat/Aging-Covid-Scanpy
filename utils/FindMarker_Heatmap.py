import scanpy as sc
import argparse
import sys
import os
import scanyuan as scy
sys.path.append("../")

from core.utils import get_rank_group_genes,subset_by_column

parser=argparse.ArgumentParser()
parser.add_argument("-o","--outdir",type=str,default="output")
parser.add_argument("-d","--data",type=str,default=None,help="data")
parser.add_argument("-s","--subset",nargs="+",type=str,default=None,help="subset")
parser.add_argument("-c","--column",type=str,default="leiden",help="data")
parser.add_argument("-g","--groupby",type=str,default="status")

args=parser.parse_args()

if not os.path.exists(args.outdir):
   os.makedirs(args.outdir)

sc.settings.figdir=args.outdir

adata=sc.read_h5ad(args.data)
#if args.subset is not None:
#    adata=subset_by_column(adata,args.subset,args.column)
#    sc.pp.normalize_total(adata, target_sum=1e4)
if args.subset is not None:
    groups=args.subset
else:
    groups="all"
sc.tl.rank_genes_groups(adata,groups=groups,groupby=args.groupby, method='t-test',n_genes=500,corr_method="bonferroni")
sc.tl.filter_rank_genes_groups(adata,groupby=args.groupby,min_fold_change=0.25)

get_rank_group_genes(adata,outdir=args.outdir)

########
markers=get_rank_group_genes(adata,outdir=args.outdir,n_top=30,back=True)
print(markers.shape)

genes=markers.names.values.tolist()
sc.pl.heatmap(adata,var_names=genes,groupby=args.groupby,
        swap_axes=True,show=False,show_gene_labels=True,
        vmin=-3, vmax=3,save=True,figsize=[10,16])

sc.pl.rank_genes_groups_heatmap(adata, n_genes=30, 
        swap_axes=True, vmin=-3, vmax=3,
        show=False,show_gene_labels=True,save="_rank_genes_groups",figsize=[10,16],dendrogram=False)

sc.pl.rank_genes_groups_heatmap(adata, n_genes=30,log=True,
        swap_axes=True,vmin=-3, vmax=3,
        show=False,show_gene_labels=True,save="_rank_genes_groups_log",figsize=[10,16],dendrogram=False)
