import scanpy as sc
import os
import numpy as np
import pandas as pd
def subset_by_column(adata,subset,col_use="orig_cluster",invert=False):
      metadata=adata.obs
      assert col_use in metadata.columns
      #sets=metadata[col_use].values.tolist()
      adata.obs[col_use].astype("category")
      #sets=[str(s) for s in sets]
      subset=[str(s) for s in subset]
      cols=np.unique(adata.obs[col_use].values.tolist())
      print(subset)
      if invert:
          subset=[s for s in cols if s not in subset]
      #    idx=[False if s in subset else True for s in sets]
      #else:
      #    idx=[True if s in subset else False for s in sets]
      #adata_subset=adata[idx,:]
      adata_subset=adata[adata.obs[col_use].isin(subset)]
      print("after subset,adata shape is [{},{}]".format(adata_subset.shape[0],adata_subset.shape[1]))
      return adata_subset

def subset_by_cell_feature(adata,cells=None,features=None):
    data=adata.copy()
    data.obs["cells"]=data.obs_names
    Features=adata.var_names
    if cells is not None:
       data=data[data.obs["cells"].isin(cells)]
    if features is not None:
       feature_idx=[True if feature in features else False for feature in Features]
       data=data[:,feature_idx]

    return data


def get_rank_group_genes(adata,outdir="./markers",ntop=None,back=False):
    assert "rank_genes_groups" in adata.uns_keys()
    groupby=adata.uns["rank_genes_groups"]["params"]["groupby"]
    G=np.unique(adata.obs[groupby].values)
    table=[]
    for g in G:
        print("Get markers from : {}".format(g))
        t=sc.get.rank_genes_groups_df(adata,group=g,pval_cutoff=0.05,log2fc_min=0.25)
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
    if back:
        return df


def AddMetaData(data,meta,inplace=False):
    if not inplace:
        adata=data.copy()
    else: 
        adata=data
    metadata=adata.obs
    assert meta.shape[0]==metadata.shape[0]
    index=metadata.index.tolist()
    DATA=meta.loc[index,:]  #  loc accept str index
    columns=DATA.columns.tolist()
    for column in columns:
        if column in metadata.columns:
            #new_column="X_"+str(column)
            continue
        else:
            new_column=str(column)
            adata.obs[new_column]=DATA[new_column].values
    if inplace:
        return 
    else:
        return adata





