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


safeBPParam <- function(nworkers) {
    if (.Platform$OS.type=="windows") {
        BiocParallel::SerialParam()
    } else {
        BiocParallel::MulticoreParam(nworkers)
    }
}
(cores <- parallel::detectCores())
(bpParam <- safeBPParam(cores))

sce<-readRDS(args$sce)
exprs="exprs"

# # KNN 的计算在 UMAP 和 TSNE 是可以复用的
library(BiocNeighbors)
x <- assay(sce, exprs_values = exprs)
x_index <- buildKmknn(x)
(knnParam <- KmknnParam(BNINDEX=x_index))


print("### Run UMAP")
p1=proc.time()
sce <- runUMAP(sce, exprs_values=exprs,
    external_neighbors=TRUE, BPPARAM=bpParam,
    approx_pow=TRUE, init="spca", n_threads=cores, pca=50, n_sgd_threads=cores)
p2=proc.time()
time_delta<-p2[3]-p1[3]
sprintf("run UMAP use %f seconds",time_delta)

print("### Run TSNE")
p1=proc.time()
sce <- runTSNE(sce, exprs_values=exprs,
    external_neighbors=TRUE, BPPARAM=bpParam,
    approx_pow=TRUE, init="spca", n_threads=cores, pca=50, n_sgd_threads=cores)
p2=proc.time()
time_delta<-p2[3]-p1[3]
sprintf("run TSNE use %f seconds",time_delta)
# BNPARAM = KmknnParam()
#sce <- runUMAP(sce, exprs_values = "exprs")
#sce <- runTSNE(sce, exprs_values = "exprs")
saveRDS(sce,file.path(args$outdir,"sce_parallel.rds"))
