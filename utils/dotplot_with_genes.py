import scanpy as sc
import argparse
import sys
import os
import scanyuan as scy
import pandas as pd

sys.path.append("../")

from core.utils import get_rank_group_genes,subset_by_column

parser=argparse.ArgumentParser()
parser.add_argument("--outdir",type=str,default="output")
parser.add_argument("--data",type=str,default=None,help="data")
parser.add_argument("--subset",nargs="+",type=str,default=None,help="subset")
parser.add_argument("--genes",type=str,default=None,help="show genes")
parser.add_argument("--column",type=str,default="leiden",help="data")
parser.add_argument("--groupby",type=str,default="status")
parser.add_argument("--level",nargs="+",type=str,default=None)

args=parser.parse_args()

if not os.path.exists(args.outdir):
   os.makedirs(args.outdir)
sc.settings.set_figure_params(dpi_save=200)
sc.settings.figdir=args.outdir

df=pd.read_csv(args.genes,sep="\t",header=None)
marker_genes=df.loc[:,0].values.tolist()
adata=sc.read_h5ad(args.data)
if args.subset is not None:
    adata=subset_by_column(adata,args.subset,args.column)

if args.level is not None:
    adata.obs[args.groupby]=adata.obs[args.groupby].cat.reorder_categories(args.level)

# "Oranges"
sc.pl.dotplot(adata,var_names=marker_genes,groupby=args.groupby,color_map="YlGnBu",dendrogram=False,
        figsize=[16,8],show=False,save=True,smallest_dot=40,standard_scale='var')
