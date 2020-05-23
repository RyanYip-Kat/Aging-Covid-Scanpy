import scanpy as sc
import argparse
import sys
import os
import scanyuan as scy
sys.path.append("../")

from core.utils import get_rank_group_genes,subset_by_column

parser=argparse.ArgumentParser()
parser.add_argument("-o","--outdir",type=str,default="output")

args=parser.parse_args()

if not os.path.exists(args.outdir):
   os.makedirs(args.outdir)

sc.settings.figdir=args.outdir

adata=sc.read_h5ad("output/aging-11/DC/cDC2-1-2/adata.h5ad")
sc.pp.normalize_total(adata, target_sum=1e4)


marker_genes=["CLEC12A", "DUSP6","TXNIP","VCAN","NEAT1","CFD","EIF3A","LAMP1","CD36", "CST3", "VIM","B2M","PSMA2","DUSP2","IFITM1","DDIT4","FCER1A","EEF1G","RGS10","HLA-DQA2"]
sc.pl.stacked_violin(adata, marker_genes, groupby="new_clusters",show=False,show_gene_labels=True,save="_cDC12_1",figsize=[16,8],dendrogram=False)
sc.pl.stacked_violin(adata, marker_genes, groupby="new_clusters",show=False,show_gene_labels=True,save="_cDC12_2",figsize=[8,16],dendrogram=False,swap_axes=True)
scy.stacked_violin_t(adata, marker_genes, figsize=[16,8], groupby='new_clusters',show=False,show_gene_labels=True,use_raw=True,save="_cDC12_3")
