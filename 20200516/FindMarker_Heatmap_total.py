import scanpy as sc
import argparse
import sys
import os
sys.path.append("../")

from core.utils import get_rank_group_genes,subset_by_column

parser=argparse.ArgumentParser()
parser.add_argument("-o","--outdir",type=str,default="output")

args=parser.parse_args()

if not os.path.exists(args.outdir):
   os.makedirs(args.outdir)

sc.settings.figdir=args.outdir

adata=sc.read_h5ad("../20200513/output/aging-11/Total/adata_doublet.h5ad")
adata.obs["celltype"]=adata.obs["celltype"].cat.reorder_categories(["TC","NK","BC","Mono","DC","RBC","MEG"])
sc.tl.rank_genes_groups(adata,groupby="celltype", method='t-test',n_genes=500,corr_method="bonferroni")
get_rank_group_genes(adata,args.outdir)

########
markers=get_rank_group_genes(adata,args.outdir,15,True)
print(markers.shape)

genes=markers.names.values.tolist()
sc.pl.heatmap(adata,var_names=genes,groupby="celltype",
        swap_axes=True,show=False,show_gene_labels=True,
        vmin=-3, vmax=3,save=True,figsize=[10,16])

sc.pl.rank_genes_groups_heatmap(adata, n_genes=15, 
        swap_axes=True, vmin=-3, vmax=3,
        show=False,show_gene_labels=True,save="_rank_genes_groups",figsize=[10,16],dendrogram=False)



