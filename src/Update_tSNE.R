library(argparse)
library(stringr)
library(flowCore)
library(FlowSOM)
library(Rtsne)
library(ggplot2)
library(sva)

parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--outdir",
                    type="character",
                    default="output",
                    help="the path to save result")


parser$add_argument("--sce",
                    type="character",
                    default="")


parser$add_argument("--tSNE",
                    type="character",
                    default="tSNE file")

args <- parser$parse_args()
if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}

sce<-readRDS(args$sce)
DATA<-readRDS(args$tSNE)

DATA<-DATA[,c(1,2)]
DATA<-DATA[colnames(sce),]

print("###  Update tSNE")
reducedDim(sce,"TSNE")=DATA
colData(sce)$metacluster=DATA[,"cluster"]

saveRDS(sce,file.path(args$outdir,"sce_update.rds"))
