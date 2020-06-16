import scanpy as sc
import harmonypy as hm

def process(adata,key=None,method="harmony",use_hv=True):
      # method : harmony,combat
      print("**** data transform ****")
      sc.pp.normalize_total(adata, target_sum=1e4)
      sc.pp.log1p(adata)
      adata.raw = adata

      sc.pp.highly_variable_genes(adata,n_top_genes=5000)
      if use_hv:
          adata = adata[:, adata.var.highly_variable]
      print("**** scale data ****")
      sc.pp.scale(adata, max_value=10)

      if key is not None:
          assert key in adata.obs.columns
          print("**** correct batch effect ****")
          #if "harmony" in adata.uns_keys() 
          #sc.pp.neighbors(adata,use_rep="X_harmony")
          #print("**** run tsne ****")
          #sc.tl.tsne(adata,use_rep="X_harmony",learning_rate=200)
          if method=="harmony":
              adata=run_harmony(adata,batch_key=key)
          elif method=="combat":
              sc.pp.combat(adata,key=key)
          else:
              raise ValueError("Method must be one of hamony or combat")

          var_genes=adata.uns["HVG"]
          print("There are {} var genes".format(len(var_genes)))
          #sc.external.pp.mnn_correct(adata,var_subset=var_genes,batch_key=key)
          #sc.external.pp.bbknn(adata,batch_key=key,approx=True)
      return adata

def reduction(adata,resolution=1.5):
      if "X_pca" not in adata.obsm_keys():
         print("Run PCA first now!!!")
         sc.tl.pca(adata,n_comps=n_pca,svd_solver='arpack')
         sc.pl.pca_variance_ratio(adata,show=False, log=True,save=True)

      use_rep="X_harmony" if "X_harmony" in adata.obsm_keys() else "X_pca"

      print("**** run neighborhood graph ****")
      sc.pp.neighbors(adata, n_neighbors=20,use_rep=use_rep)

      print("**** run tsne ****")
      sc.tl.tsne(adata,use_rep=use_rep,learning_rate=200)

      print("**** clustering ****")
      sc.tl.leiden(adata,resolution=resolution)
      #sc.tl.louvain(adata,resolution=resolution)

      print("**** run umap ****")
      sc.tl.paga(adata)  # remove `plot=False` if you want to see the coarse-grained graph
      sc.pl.paga(adata, plot=False)
      sc.tl.umap(adata, init_pos='paga')

      return adata



def run_harmony(adata,n_pca=30,batch_key="status",max_iter=10):
    data=adata.copy()
    #if "X_pca" not in data.obsm_keys():
    #    print("Run PCA first now!!!")
    #    sc.tl.pca(data,n_comps=n_pca)
    sc.tl.pca(data,n_comps=n_pca)
    X=data.obsm["X_pca"]
    if n_pca > X.shape[1]:
        n_pca=X.shape[1]
    X=X[:,:n_pca]
    meta_data=data.obs
    assert batch_key in meta_data.columns
    vars_use=[batch_key]
    print("*** run harmony ***")
    ho = hm.run_harmony(X, meta_data, vars_use,max_iter_kmeans=max_iter)
    X_harmony=ho.Z_corr.transpose()
    
    print("*** add harmony ***")
    data.uns["harmony"]={}
    data.uns["harmony"]["params"]=ho.__dict__
    data.obsm["X_harmony"]=X_harmony
    return data

