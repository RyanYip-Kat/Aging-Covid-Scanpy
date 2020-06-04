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
parser.add_argument("-n","--n_genes",type=int,default=20,help="data")
parser.add_argument("-s","--subset",nargs="+",type=str,default=None,help="subset")
parser.add_argument("-c","--column",type=str,default="leiden",help="data")
parser.add_argument("-g","--groupby",type=str,default="status")

args=parser.parse_args()

if not os.path.exists(args.outdir):
   os.makedirs(args.outdir)
sc.settings.set_figure_params(dpi_save=200)
sc.settings.figdir=args.outdir

adata=sc.read_h5ad(args.data)
if args.subset is not None:
    adata=subset_by_column(adata,args.subset,args.column)
    adata=adata.raw.to_adata()
    adata.raw=adata
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.highly_variable_genes(adata,n_top_genes=5000)
    adata= adata[:, adata.var.highly_variable]
    sc.pp.scale(adata, max_value=10)
    adata.obs[args.column]=adata.obs[args.column].astype("category")
    adata.obs[args.column]=adata.obs[args.column].cat.reorder_categories(args.subset)
#if args.subset is not None:
#    groups=args.subset
#else:
#    groups="all"
#sc.tl.rank_genes_groups(adata,groups=groups,groupby=args.groupby, method='t-test',n_genes=500,corr_method="bonferroni")
sc.tl.rank_genes_groups(adata,groupby=args.groupby, method='t-test',n_genes=500,corr_method="bonferroni")
get_rank_group_genes(adata,pval=None,fc=0.25,outdir=args.outdir)

########

sc.pl.rank_genes_groups_heatmap(adata,n_genes=args.n_genes, 
        swap_axes=True, vmin=-3, vmax=3,
        show=False,show_gene_labels=True,save="_rank_genes_groups_1",figsize=[24,24],dendrogram=False)

sc.pl.rank_genes_groups_heatmap(adata,n_genes=args.n_genes,
        swap_axes=False, vmin=-3, vmax=3,
        show=False,show_gene_labels=True,save="_rank_genes_groups_2",figsize=[24,24],dendrogram=False)

sc.pl.rank_genes_groups_heatmap(adata,n_genes=args.n_genes,log=True,
        swap_axes=True,vmin=-3, vmax=3,
        show=False,show_gene_labels=True,save="_rank_genes_groups_log",figsize=[24,24],dendrogram=False)
