import stream as st
import scanpy as sc
import matplotlib.pyplot as plt

import numpy  as np
import pandas as pd
import os
import argparse

from src.core import *

#####################
parser=argparse.ArgumentParser()
parser.add_argument("-o","--outdir",type=str,default=None)
parser.add_argument("-w","--workdir",type=str,default=None)
parser.add_argument("-d","--data",type=str,default=None)
parser.add_argument("-n","--n_jobs",type=int,default=8)
args=parser.parse_args()
#####################
sc.logging.print_versions()

####################
workdir=args.workdir
files=os.listdir(workdir)
if "adata.h5ad" not in files \
        or "cell_label.tsv" not in files \
        or "cell_label_color.tsv" not in files:
            raise ValueError("Check these files in workdir")

if args.outdir is None:
    outdir=workdir
else:
    outdir=args.outdir
    if not os.path.exists(outdir):
        os.makedirs(outdir)
##################
data=os.path.join(workdir,"adata.h5ad")
cell_label=os.path.join(workdir,"cell_label.tsv")
cell_label_color=os.path.join(workdir,"cell_label_color.tsv")

###################
print("### Loading data")
adata=read(data,file_format="h5ad",workdir=outdir)
adata=adata.raw.to_adata() # use raw data
adata.obs.head()
add_cell_labels(adata,file_name=cell_label)
add_cell_colors(adata,file_name=cell_label_color)
###################

filter_genes(adata,min_num_cells = 5)
print("### Select top pcs")
select_top_principal_components(adata,n_pc=30,first_pc=True,save_fig=True,fig_size=(12,12))

print("### Select variable genes")
select_variable_genes(adata,loess_frac=0.01, percentile=95,n_jobs=args.n_jobs,save_fig=True,fig_size=(12,12))

print("### dimension reduction")
dimension_reduction(adata,n_neighbors=20,n_components=2,n_jobs=args.n_jobs,method="se",feature="var_genes") # top_pcs

st.plot_dimension_reduction(adata,save_fig=True,fig_size=(12,12))
print("### elastic_principal_graph")
st.seed_elastic_principal_graph(adata,n_clusters=10)

st.plot_branches(adata,save_fig=True,fig_size=(12,12))
st.plot_branches_with_cells(adata,save_fig=True,fig_size=(12,12),fig_legend=False)

print("### elastic_principal_graph")
st.elastic_principal_graph(adata,epg_alpha=0.02,epg_mu=0.1,epg_lambda=0.02)

st.plot_branches(adata,save_fig=True,fig_size=(12,12))
st.plot_branches_with_cells(adata,save_fig=True,fig_size=(12,12))

print("### optimize_branching")
st.optimize_branching(adata,epg_alpha=0.02,epg_mu=0.1,epg_lambda=0.01)
st.plot_branches(adata,save_fig=True,fig_size=(12,12))
st.plot_branches_with_cells(adata,save_fig=True,fig_size=(12,12))


print("###Extend leaf branch to reach further cells")
st.extend_elastic_principal_graph(adata)
st.plot_branches(adata,save_fig=True,fig_size=(12,12))
st.plot_branches_with_cells(adata,save_fig=True,fig_size=(12,12))

st.plot_flat_tree(adata,fig_legend=False,save_fig=True,fig_size=(12,12))

print("### Detect genes")
try:
    detect_transistion_genes(adata,root='S0',n_jobs=args.n_jobs)

except:
    print("detect_transistion_genes failed!!!")

try:
    #from src.core import detect_de_genes,detect_leaf_genes
    detect_leaf_genes(adata,root='S0',cutoff_zscore=1.0,n_jobs=args.n_jobs)
except:
    print("detect_leaf_genes failed!!!")

try:
    detect_de_genes(adata,root='S0',n_jobs=args.n_jobs)
except:
    print("detect_de_genes failed!!!")

st.write(adata,file_name="adata_stream.pkl")
