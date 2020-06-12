import scanpy as sc
import argparse
import sys
import os
import seaborn as sns
sys.path.append("../")
from core.utils import subset_by_column,subset_by_cell_feature
from collections import Counter
from core.plotting import violin_hue,violin_sns

parser=argparse.ArgumentParser()
parser.add_argument("-o","--outdir",type=str,default="output")
parser.add_argument("-d","--data",type=str,default=None,help="data")
parser.add_argument("-g","--genes",nargs="+",type=str,default=None)
parser.add_argument("-r","--orient",type=str,default="h",help="")
parser.add_argument("-i","--invert", action='store_true', default=False)
parser.add_argument("-s","--subset",nargs="+",type=str,default=None)
parser.add_argument("-c","--column",type=str,default="leiden")
parser.add_argument("--groupby",type=str,default="status",help="")
args=parser.parse_args()

if not os.path.exists(args.outdir):
   os.makedirs(args.outdir)

genes=args.genes
adata=sc.read_h5ad(args.data)

if args.subset is not None:
    print("*** subset")
    print("Subset data by {}".format(args.column))
    adata=subset_by_column(adata,args.subset,args.column,invert=args.invert)

if args.groupby=="status":
    adata.obs[args.groupby]=adata.obs[args.groupby].cat.reorder_categories(["YA","AA"])

for gene in genes:
    #sc.pl.violin(adata,gene,groupby="status",show=False,save="_"+gene+"_status.pdf")
    violin_sns(adata,gene_name=gene,save=True,groupby=args.groupby,figdir=args.outdir,orient=args.orient,strip=True)
