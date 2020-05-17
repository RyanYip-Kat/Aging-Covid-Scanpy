import scanpy as sc

def process(adata,key=None,use_hv=True):
      print("**** data transform ****")
      sc.pp.normalize_total(adata, target_sum=1e4)
      sc.pp.log1p(adata)
      adata.raw = adata

      sc.pp.highly_variable_genes(adata,n_top_genes=5000)
      if use_hv:
          adata = adata[:, adata.var.highly_variable]
      print("**** scale data ****")
      sc.pp.scale(adata, max_value=10)

      #sc.pp.combat(adata,key="idents")
      if key is not None:
          assert key in adata.obs.columns
          print("**** correct batch effect ****")
          sc.pp.combat(adata,key=key)
          var_genes=adata.uns["HVG"]
          print("There are {} var genes".format(len(var_genes)))
          #sc.external.pp.mnn_correct(adata,var_subset=var_genes,batch_key=key)
          #sc.external.pp.bbknn(adata,batch_key=key,approx=True)
      print(adata)
      return adata

def reduction(adata,resolution=1.5):
      print("*** Principal component analysis ***")
      print("**** run pca ****")
      sc.tl.pca(adata, svd_solver='arpack')
      sc.pl.pca_variance_ratio(adata, log=True,save=True)
      print("**** run neighborhood graph ****")
      sc.pp.neighbors(adata, n_neighbors=20, n_pcs=30)

      print("**** run tsne ****")
      sc.tl.tsne(adata,n_pcs=30,learning_rate=200)

      print("**** clustering ****")
      sc.tl.leiden(adata,resolution=resolution)
      sc.tl.louvain(adata,resolution=resolution)

      print("**** run umap ****")
      sc.tl.paga(adata)  # remove `plot=False` if you want to see the coarse-grained graph
      sc.pl.paga(adata, plot=False)
      sc.tl.umap(adata, init_pos='paga')

      return adata
