import os
import scanpy as sc
import scirpy as ir
import argparse

import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt

from collections import Counter

class Model(object):
   def __init__(self,path,tcr_path=None,n_top_genes=5000,outdir="./scanpy_result"):
      self._path=path
      self._outdir=outdir
      self._tcr_path=tcr_path

      self._n_top_genes=n_top_genes

      self.graphclust_path=os.path.join(self._path,"outs/analysis/clustering/graphclust/clusters.csv")
      self.km_path=os.path.join(self._path,"analysis/clustering/kmeans_10_clusters/clusters.csv")
      self.mtx_path=os.path.join(self._path,"outs/filtered_feature_bc_matrix")

   def _load_data(self):
      adata=sc.read_10x_mtx(self.mtx_path,cache=True)
      orig_cluster=pd.read_csv(self.graphclust_path)
      km_cluster=pd.read_csv(self.km_path)

      print("*** Add original message")
      adata.obs["orig_cluster"]=orig_cluster.Cluster.to_list()
      adata.obs["km_cluster"]=km_cluster.Cluster.to_list()
       
      cells=adata.obs_names.to_list()
      idents=[cell.split("-")[1] for cell in cells]
      adata.obs["idents"]=idents

      self._n_cells=adata.shape[0]
      self._n_features=adata.shape[1]
      print("origin adata shape is [{},{}]".format(self._n_cells,self._n_features))
      adata=self._remove(adata)  # remove dropout features

      if self._tcr_path is not None:
            if os.path.exists(self._tcr_path):
               adata.uns["tcr_path"]=self._tcr_path
               try:
                  adata_tcr=ir.read_10x_vdj(self._tcr_path)
                  print("*** Merge adata with tcr")
                  ir.pp.merge_with_tcr(adata,adata_tcr)

                  print("*** chain pairing")
                  ir.tl.chain_pairing(adata)

                  ir.pl.group_abundance(adata, groupby="chain_pairing", target_col="orig_cluster")
                  plt.savefig(os.path.join(args.outdir,"chain_pairing.pdf"))
                  plt.close()

                  print("Fraction of cells with more than one pair of TCRs: {:.2f}".format(np.sum(
                     adata.obs["chain_pairing"].isin(["Extra beta", "Extra alpha", "Two full chains"])) / adata.n_obs))
               except:
                  print("Add TCR Invalid!!!")

      #if self._subset is not None:
      #  adata=self._select(adata)
      #self.adata=adata
      adata.uns["path"]=self._path
      return adata


   
