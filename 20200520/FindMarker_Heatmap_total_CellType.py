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

adata=sc.read_h5ad("../20200517/output/aging-11/Total/CellType/adata_CellType.h5ad")
levels=["CD4 Naive","CD4 Tcm","CD4 Tem","CD4 Treg","CD4 Tex",
        "CD8 Naive","CD8 Tem","CD8 CTL","CD8 Tex",
        "CD4-CD8-","CD4+CD8+","T-mito","NK1","NK2","NK3",
        "Naive","Memory","ASC","ABC","CD14","CD16","Intermed","pDC","cDC1","cDC2"]



print(Counter(adata.obs["CellType"]))
adata.obs["CellType"]=adata.obs["CellType"].cat.reorder_categories(levels)
adata=adata[adata.obs["CellType"].isin(["CD4 Naive","Naive","CD14","cDC2"])]
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.scale(adata, max_value=10)

sc.tl.rank_genes_groups(adata,groupby="CellType",n_genes=500,corr_method="bonferroni")
sc.tl.filter_rank_genes_groups(adata,groupby="CellType",min_fold_change=0.25)
get_rank_group_genes(adata,args.outdir)

########
markers=get_rank_group_genes(adata,args.outdir,20,True)
print(markers.shape)

sc.pl.rank_genes_groups_heatmap(adata, n_genes=25, 
        swap_axes=True, vmin=-3, vmax=3,
        show=False,show_gene_labels=True,save="_rank_genes_groups_1",figsize=[16,20],dendrogram=False)

sc.pl.rank_genes_groups_heatmap(adata, n_genes=25,
        vmin=-3, vmax=3,
        show=False,show_gene_labels=True,save="_rank_genes_groups_2",figsize=[20,16],dendrogram=False)

sc.pl.rank_genes_groups_heatmap(adata, n_genes=25,key="rank_genes_groups_filtered",
        swap_axes=True, vmin=-3, vmax=3,
        show=False,show_gene_labels=True,save="_rank_genes_groups_filtered_1",figsize=[16,20],dendrogram=False)

sc.pl.rank_genes_groups_heatmap(adata, n_genes=25,key="rank_genes_groups_filtered",
        vmin=-3, vmax=3,
        show=False,show_gene_labels=True,save="_rank_genes_groups_filtered_2",figsize=[20,16],dendrogram=False)


sc.pl.rank_genes_groups_dotplot(adata, n_genes=25,show=False,save="_rank_genes_groups_1",dendrogram=False,figsize=[20,16])



