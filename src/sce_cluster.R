library(CATALYST)
library(cowplot)
library(scater)
library(SingleCellExperiment)
library(argparse)
library(stringr)

parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--outdir",
                    type="character",
                    default="output",
                    help="the path to save result")


parser$add_argument("--sce",
		    type="character",
		    default="")

args <- parser$parse_args()
if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}

sce<-readRDS(args$sce)
safeBPParam <- function(nworkers) {
    if (.Platform$OS.type=="windows") {
        BiocParallel::SerialParam()
    } else {
        BiocParallel::MulticoreParam(nworkers)
    }
}
(cores <- parallel::detectCores())
(bpParam <- safeBPParam(cores))

lineage_markers <- rownames(sce) 
rownames(sce)=str_replace(lineage_markers,"-","_")
rowData(sce)$antigen=str_replace(lineage_markers,"-","_")
#lineage_markers<-readRDS(args$path,"antigen.rds")
sce <- cluster(sce, features = lineage_markers, 
    xdim = 8, ydim = 8, maxK = 50, verbose = TRUE, seed = 1)     

saveRDS(sce,file.path(args$outdir,"sce.rds"))
