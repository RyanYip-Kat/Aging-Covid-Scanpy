import scanpy as sc
import pandas as pd
import argparse
import sys
import os
sys.path.append("../")

from core.utils import get_rank_group_genes,subset_by_column
from collections import Counter
parser=argparse.ArgumentParser()
parser.add_argument("-o","--outdir",type=str,default="output")

args=parser.parse_args()

if not os.path.exists(args.outdir):
   os.makedirs(args.outdir)

sc.settings.figdir=args.outdir

adata=sc.read_h5ad("../20200517/output/aging-11/Total/CellType/adata_CellType.h5ad")
levels=["CD4 Naive","CD4 Tcm","CD4 Tem","CD4 Treg","CD4 Tex",
        "CD8 Naive","CD8 Tem","CD8 CTL","CD8 Tex",
        "CD4-CD8-","CD4+CD8+","T-mito","NK1","NK2","NK3",
        "Naive","Memory","ASC","ABC","CD14","CD16","Intermed","pDC","cDC1","cDC2"]



print(Counter(adata.obs["CellType"]))
adata.obs["CellType"]=adata.obs["CellType"].cat.reorder_categories(levels)
markers=pd.read_csv("../20200517/output/aging-11/Total/CellType/top4_markers.csv")
#markers["cluster"]=markers["cluster"].astype("category")
#markers["cluster"]=markers["cluster"].cat.reorder_categories(levels)
#markers=markers.sort_values(by=["cluster"])
#genes=markers["names"].values.tolist()
genes=pd.read_csv("genes.txt",sep="\t",header=None)
genes=genes.loc[:,0].values.tolist()

sc.pl.heatmap(adata,use_raw=True,var_names=genes,groupby="CellType",show_gene_labels=True,
        show=False,swap_axes=True, vmin=-3, vmax=3,
        save="genes_groups_1",figsize=[16,20],dendrogram=False)

sc.pl.heatmap(adata,use_raw=True,var_names=genes,groupby="CellType",show_gene_labels=True,
        show=False,swap_axes=False, vmin=-3, vmax=3,
        save="genes_groups_2",figsize=[20,16],dendrogram=False)






