import scanpy as sc
import os
import numpy as np
import pandas as pd
def subset_by_column(adata,subset,col_use="orig_cluster",invert=False):
      metadata=adata.obs
      assert col_use in metadata.columns
      sets=metadata[col_use].values.tolist()
      sets=[str(s) for s in sets]
      subset=[str(s) for s in subset]
      print(subset)
      if invert:
          idx=[False if s in subset else True for s in sets]
      else:
          idx=[True if s in subset else False for s in sets]
      adata_subset=adata[idx,:]
      print("after subset,adata shape is [{},{}]".format(adata_subset.shape[0],adata_subset.shape[1]))
      return adata_subset

def subset_by_cell_feature(adata,cells=None,features=None):
    Cells=adata.obs_names
    Features=adata.var_names
    if cells is not None:
       cell_idx=[True if cell in cells else False for cell in Cells]
       adata=adata[cell_idx,:]
    if features is not None:
       feature_idx=[True if feature in features else False for feature in Features]
       adata=adata[:,feature_idx]

    return adata


def get_rank_group_genes(adata,outdir="./markers",ntop=None):
    assert "rank_genes_groups" in adata.uns_keys()
    groupby=adata.uns["rank_genes_groups"]["params"]["groupby"]
    G=np.unique(adata.obs[groupby].values)
    table=[]
    for g in G:
        print("Get markers from : {}".format(g))
        t=sc.get.rank_genes_groups_df(adata,group=g)
        t["cluster"]=g
        t=t.sort_values("logfoldchanges",ascending=False)
        if ntop is not None:
            t=t[:ntop]

        table.append(t)
     
    df=pd.concat(table,axis=0)
    if ntop is not None:
        filename=os.path.join(outdir,str(ntop)+"_markers.csv")
    else:
        filename=os.path.join(outdir,"markers.csv")
    print("Save markers in {}".format(filename))
    df.to_csv(filename,sep=",",index=False)
