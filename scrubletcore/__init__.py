# some technical stuff
import sys
# the actual API
from . import tools as tl
from . import preprocessing as pp
from . import plotting as pl
from . import datasets, logging, queries, external

from anndata import AnnData
from anndata import read_h5ad, read_csv, read_excel, read_hdf, read_loom, read_mtx, read_text, read_umi_tools
from .readwrite import read, read_10x_h5, read_10x_mtx, write
from .neighbors import Neighbors
from ._settings import settings

set_figure_params = settings.set_figure_params

# has to be done at the end, after everything has been imported
