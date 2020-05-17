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
parser.add_argument("-g","--groupby",type=str,default="status")

args=parser.parse_args()

if not os.path.exists(args.outdir):
   os.makedirs(args.outdir)

sc.settings.figdir=args.outdir

adata=sc.read_h5ad(args.data)
if args.subset is not None:
    adata=subset_by_column(adata,args.subset,args.column)
sc.tl.rank_genes_groups(adata,groupby=args.groupby, method='t-test',n_genes=500,corr_method="bonferroni")
get_rank_group_genes(adata,args.outdir)

########
markers=get_rank_group_genes(adata,args.outdir,15,True)
print(markers.shape)

genes=markers.names.values.tolist()
sc.pl.heatmap(adata,var_names=genes,groupby=args.groupby,
        swap_axes=True,show=False,show_gene_labels=True,
        vmin=-3, vmax=3,save=True,figsize=[10,16])

sc.pl.rank_genes_groups_heatmap(adata, n_genes=15, 
        swap_axes=True, vmin=-3, vmax=3,
        show=False,show_gene_labels=True,save="_rank_genes_groups",figsize=[10,16],dendrogram=False)



