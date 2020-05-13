
import stream as st
import scanpy as sc
import matplotlib.pyplot as plt

import numpy  as np
import pandas as pd
import os

from src.core import *


st.__version__


# ### Loading data from scanpy object


adata=st.read(file_name="stream_result/BC/adata.h5ad",workdir="./stream_result",file_format="h5ad")


# ### Add meta data

# In[5]:


st.add_cell_labels(adata,file_name="stream_result/BC/cell_label.tsv")
st.add_cell_colors(adata,file_name="stream_result/BC/cell_label_color.tsv")


# In[6]:


adata.uns["label_color"]


# In[7]:


st.select_top_principal_components(adata,n_pc=20,first_pc=True)


# In[8]:


st.select_variable_genes(adata,n_genes=500)


# ### Run Dimension reduction

# In[9]:


st.dimension_reduction(adata,n_neighbors=20,n_components=2,n_jobs=8,method="se",feature="top_pcs")


# In[10]:


st.plot_dimension_reduction(adata)
st.plot_visualization_2D(adata,use_precomputed=False)


# In[20]:


st.seed_elastic_principal_graph(adata,n_clusters=15)


# In[21]:


st.plot_branches(adata)
st.plot_branches_with_cells(adata,fig_legend=False)


# In[14]:


#st.seed_elastic_principal_graph(adata,n_clusters=10)
st.elastic_principal_graph(adata,epg_alpha=0.02,epg_mu=0.1,epg_lambda=0.02,epg_n_processes=8)


# ### optional step

# In[16]:


### optional step
#st.optimize_branching(adata,epg_alpha=0.02,epg_mu=0.1,epg_lambda=0.01)
#st.plot_branches(adata)
#st.plot_branches_with_cells(adata)


# In[22]:


st.subwaymap_plot(adata,root='S0',fig_legend_ncol=5) 


# In[23]:


st.stream_plot(adata,root='S0',fig_legend_ncol=5,fig_size=(12,12))

