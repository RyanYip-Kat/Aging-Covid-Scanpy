import argparse
import os
import stream as st
import scanpy as sc
import numpy as np

parser=argparse.ArgumentParser()
parser.add_argument("-o","--outdir",type=str,default=None)

args=parser.parse_args()

adata=st.read("ouput/aging-xxx/subset/adata.h5ad",file_format="h5ad",workdir=os.path.join(args.outdir,"stream_result"))
adata.obs["label"]=adata.obs["louvain"]

label_color={"0": "#FF0000",
             "1": "#836FFF",
             "2": "#0000FF",
             "3": "#C6E2FF",
             "4": "#548B54",
             "5": "#00FF00",
             "6": "#FFF68F",
             "7": "#8B864E",
             "8": "#FFFF00",
             "9": "#FFD700",
             "10": "#8B658B",
             "11": "#FF6A6A",
             "12": "#FFD39B",
             "13": "#EE2C2C",
             "14": "#BF3EFF"
             }
             
adata.obs["label_color"] = [label_color[x] for x in adata.obs["label"]]
# st.select_variable_genes(adata,loess_frac=0.01, percentile=90)

print("### dimension_reduction")
st.dimension_reduction(adata,n_neighbors=50,feature="var_genes",method="se",n_components=3,n_jobs=8)
st.plot_dimension_reduction(adata,save_fig=True,fig_size=[16,16])

print("### seed_elastic_principal_graph")
st.seed_elastic_principal_graph(adata,n_clusters=10)
st.plot_branches(adata,save_fig=True,fig_size=[16,16],fig_name="branches_seed.pdf")
st.plot_branches_with_cells(adata,fig_legend=False,save_fig=True,fig_size=[16,16],fig_name="branches_with_cells_seed.pdf")

print("### elastic_principal_graph")
st.elastic_principal_graph(adata,epg_alpha=0.02,epg_mu=0.1,epg_lambda=0.02)
st.plot_branches(adata,save_fig=True,fig_size=[16,16],fig_name="branches_elastic.pdf")
st.plot_branches_with_cells(adata,fig_legend=False,save_fig=True,fig_size=[16,16],fig_name="branches_with_cells_elastic.pdf")

###Extend leaf branch to reach further cells 
print("### Extend leaf branch to reach further cells ")
st.extend_elastic_principal_graph(adata)
st.plot_branches(adata,save_fig=True,fig_size=[16,16],fig_name="branches_extend_elastic.pdf")
st.plot_branches_with_cells(adata,fig_legend=False,save_fig=True,fig_size=[16,16],fig_name="branches_with_cells_extend_elastic.pdf")

print("### plot_flat_tree")
st.plot_flat_tree(adata,fig_legend=False,save_fig=True,fig_size=[16,16])

print("### write tsv file for stream command line")
X=adata.X
X.tofile(os.path.join(args.outdir,"adata.tsv"),sep="\t")
cell_label=adata.obs["label"]
cell_label.to_csv("stream_result/cell_label.tsv",sep="\t",index=False)

cell_label_color=adata.obs["label_color"]
cell_label_color.to_csv("stream_result/cell_label_color.tsv",sep="\t",index=False)

adata.write(os.path.join(args.outdir,"adata_stream.h5ad"), compression='gzip')



