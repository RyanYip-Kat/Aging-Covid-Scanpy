import scanpy as sc
import argparse
import sys
import os
import seaborn as sns
sys.path.append("../")

from core.plotting import violin_hue,violin_sns

parser=argparse.ArgumentParser()
parser.add_argument("-o","--outdir",type=str,default="output")
parser.add_argument("-d","--data",type=str,default=None,help="data")
parser.add_argument("-g","--genes",nargs="+",type=str,default=None)
parser.add_argument("-c","--orient",type=str,default="h",help="")
args=parser.parse_args()

if not os.path.exists(args.outdir):
   os.makedirs(args.outdir)

genes=args.genes
adata=sc.read_h5ad(args.data)

for gene in genes:
    #sc.pl.violin(adata,gene,groupby="status",show=False,save="_"+gene+"_status.pdf")
    violin_sns(adata,gene_name=gene,save=True,figdir=args.outdir,orient=args.orient,strip=True)
