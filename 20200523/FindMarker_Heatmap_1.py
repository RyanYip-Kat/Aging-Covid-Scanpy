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


marker_genes_1=["CLEC12A","TXNIP","DUSP6", "VCAN","NEAT1", "CFD", "EIF3A", "LAMP1", "CD36", "MS4A7"]
marker_genes_2=["B2M","PSMA2","DUSP2","IFITM1","DDIT4","FCER1A","EEF1G","RGS10","HLA-DQA2","IL1B"]
sc.pl.stacked_violin(adata, marker_genes_1, groupby="new_clusters",show=False,show_gene_labels=True,save="_cDC12_1",figsize=[8,16],dendrogram=False,swap_axes=True)
sc.pl.stacked_violin(adata, marker_genes_2, groupby="new_clusters",show=False,show_gene_labels=True,save="_cDC12_2",figsize=[8,16],dendrogram=False,swap_axes=True)
