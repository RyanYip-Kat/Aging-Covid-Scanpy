import scanpy as sc
import os
import sys
import argparse
from collections import Counter

sys.path.append("../")
from core.core import Model
from core.ScanpyConfig import *
from core.Config import get_args
from core.scrubletDoublet import DoubletModule

if __name__=="__main__":
    args=get_args()
    sc.settings.n_jobs=args.n_jobs

    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)
    # path : /home/ye/Work/BioAligment/10X/output/aging-11-YA
    # tcr_path : /home/ye/Work/BioAligment/10X/output/aging-11-YA-tCR
    net=Model(args.path,tcr_path=args.tcr,use_harmony=True,n_top_genes=2000)
    print(net)
    adata=net._load_data()
    print(adata.shape)
    print(adata.shape)

    metadata=adata.obs
    idents=metadata.idents
    adata.obs["idents"]=idents

    adata=net._process(adata,batch_key="idents",remove_doublet=True)
    #adata=net._process(adata)
    adata=net._reduction(adata)
    adata=net._FindMarkers(adata)
    
    print("### Save")
    adata.write(os.path.join(args.outdir,'adata.h5ad'), compression='gzip')

