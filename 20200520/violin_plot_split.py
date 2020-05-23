import scanpy as sc
import os
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.sparse import issparse
from matplotlib import pyplot as pl
from anndata import AnnData
from typing import Union, Optional, List, Sequence, Iterable
from matplotlib.colors import Colormap


def violin_hue(
    adata: AnnData,
    groupby: Optional[Sequence[str]] = None,
    gene_name: Optional[Iterable[str]] = None,
    use_raw: Optional[bool] = None,
    hue : Optional[str] = None,
    hue_order: Optional[str] = None,
    figdir:Optional[str] ="./figdir",
    split: bool = True,
    scale: str = 'width',
    strip: bool = True,
    jitter: Union[int, float, bool] = True,
    size: int = 1,
    show: Optional[bool] = None,
    save: Optional[bool] = None,
):
    if use_raw is None:
        use_raw = bool(adata.uns[key]['params']['use_raw'])
    if adata.raw is not None and use_raw:
        X= adata.raw[:, gene_name].X
                
    else:
        X= adata[:, gene_name].X
    if issparse(X):
        X = X.toarray()
    df = pd.DataFrame(X, index=adata.obs_names, columns=["expression"])
    assert groupby in adata.obs.columns
    
    df["group"]=adata.obs[groupby]
    assert hue in adata.obs.columns
    df["hue"]=adata.obs[hue]
    df["hue"]=df["hue"].cat.reorder_categories(hue_order)

    _ax = sns.violinplot(x="group", y="expression", data=df, inner=None,
                             hue_order=None, hue='hue', split=split,
                             scale=scale, orient='vertical')
    if strip:
        _ax = sns.stripplot(x="group", y="expression", data=df,
                                hue='hue', dodge=True,hue_order=None,
                                jitter=jitter, color='black', size=size,ax=_ax)
    _ax.set_xlabel('group')
    _ax.set_ylabel('expression')
    _ax.set_title(gene_name)
    _ax.legend_.remove()
    if show:
        pl.show()
    if save:
        if not os.path.exists(figdir):
            os.makedirs(figdir)
        filename=os.path.join(figdir,gene_name+"_violin_hue.pdf")
        pl.savefig(filename, dpi=80, bbox_inches='tight')
