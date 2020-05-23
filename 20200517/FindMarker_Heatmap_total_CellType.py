import scanpy as sc
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

adata=sc.read_h5ad("output/aging-11/Total/CellType/adata_CellType.h5ad")
levels=["CD4 Naive","CD4 Tcm","CD4 Tem","CD4 Treg","CD4 Tex",
        "CD8 Naive","CD8 Tem","CD8 CTL","CD8 Tex",
        "CD4-CD8-","CD4+CD8+","T-mito","NK1","NK2","NK3",
        "Naive","Memory","ASC","ABC","CD14","CD16","Intermed","pDC","cDC1","cDC2"]



print(Counter(adata.obs["CellType"]))
adata.obs["CellType"]=adata.obs["CellType"].cat.reorder_categories(levels)
sc.tl.rank_genes_groups(adata,groupby="CellType", method='t-test',n_genes=500,corr_method="bonferroni")
get_rank_group_genes(adata,args.outdir)

########
markers=get_rank_group_genes(adata,args.outdir,8,True)
print(markers.shape)

genes=markers.names.values.tolist()
sc.pl.heatmap(adata,var_names=genes,groupby="CellType",
        swap_axes=True,show=False,show_gene_labels=True,
        vmin=-3, vmax=3,save=True,figsize=[12,16])

sc.pl.rank_genes_groups_heatmap(adata, n_genes=5, 
        swap_axes=True, vmin=-3, vmax=3,
        show=False,show_gene_labels=True,save="_rank_genes_groups_1",figsize=[16,20],dendrogram=False)

sc.pl.rank_genes_groups_heatmap(adata, n_genes=5,
        vmin=-3, vmax=3,
        show=False,show_gene_labels=True,save="_rank_genes_groups_2",figsize=[20,16],dendrogram=False)



