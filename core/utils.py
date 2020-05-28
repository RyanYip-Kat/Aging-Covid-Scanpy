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
          #subset=[s for s in cols if s not in subset]
          adata_subset=adata[~adata.obs[col_use].isin(subset)]
      else:
          adata_subset=adata[adata.obs[col_use].isin(subset)]
      print("after subset,adata shape is [{},{}]".format(adata_subset.shape[0],adata_subset.shape[1]))
      return adata_subset

def subset_by_cell_feature(adata,cells=None,features=None,invert=False):
    data=adata.copy()
    data.obs["cells"]=data.obs_names
    Features=adata.var_names
    if cells is not None:
        if invert:
            data=data[~data.obs["cells"].isin(cells)]
        else:
            data=data[data.obs["cells"].isin(cells)]
    if features is not None:
       if invert:
           feature_idx=[False if feature in features else True for feature in Features]
       else:
           feature_idx=[True if feature in features else False for feature in Features]
       data=data[:,feature_idx]

    return data


def get_rank_group_genes(adata,pval=None,fc=0.25,outdir="./markers",
        n_top=None,back=False,
        key="rank_genes_groups",prefix=""):
    # key : rank_genes_groups or rank_genes_groups_filtered
    assert "rank_genes_groups" in adata.uns_keys()
    groupby=adata.uns[key]["params"]["groupby"]
    G=np.unique(adata.obs[groupby].values)
    table=[]
    for g in G:
        try:
           print("Get markers from : {}".format(g))
           t=sc.get.rank_genes_groups_df(adata,group=g,pval_cutoff=pval,log2fc_min=fc)
           t["cluster"]=g
           t=t.sort_values("logfoldchanges",ascending=False)
           if n_top is not None:
              t=t[:n_top]

           table.append(t)
        except:
            print("group : {} has no rank".format(g))
     
    df=pd.concat(table,axis=0)
    if n_top is not None:
        filename=os.path.join(outdir,prefix+"_"+str(n_top)+"_markers.csv")
    else:
        filename=os.path.join(outdir,prefix+"_markers.csv")
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





