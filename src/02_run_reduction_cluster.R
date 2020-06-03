library(CATALYST)
library(cowplot)
library(flowCore)
library(scater)
library(SingleCellExperiment)
library(argparse)

parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--outdir",
                    type="character",
                    default="output",
                    help="the path to save result")


parser$add_argument("--path",
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


files=file.path(args$path,list.files(args$path,pattern=".fcs"))
pbmc_fs <- read.flowSet(files, transformation = FALSE, truncate_max_range = FALSE)

pbmc_panel<-readRDS(file.path(args$path,"pbmc_panel.rds"))
pbmc_md<-readRDS(file.path(args$path,"pbmc_md.rds"))
#pbmc_md$file_name=files
(sce <- prepData(pbmc_fs, pbmc_panel,pbmc_md))



print("### ExprHeatmap")
pdf(file.path(args$outdir,"plotExprHeatmap.pdf"),width=12,height=8)
plotExprHeatmap(sce, bin_anno = FALSE, row_anno = TRUE)
dev.off()

print("### Counts")
pdf(file.path(args$outdir,"Counts.pdf"),width=12,height=8)
plotCounts(sce, color_by = "condition")
dev.off()


# specify markers to use for clustering

#lineage_markers <- rownames(sce) 
lineage_markers<-readRDS(args$path,"antigen.rds")
sce <- cluster(sce, features = lineage_markers, 
    xdim = 10, ydim = 10, maxK = 20, verbose = FALSE, seed = 1)     

asinh_scale=5
DATA=assay(sce,"exprs")
DATA=asinh(DATA/asinh_scale)
assay(sce,"asinh")=DATA

exprs="asinh"

# # KNN 的计算在 UMAP 和 TSNE 是可以复用的
library(BiocNeighbors)
x <- assay(sce, exprs_values = exprs)
x_index <- buildKmknn(x)
(knnParam <- KmknnParam(BNINDEX=x_index))


print("### Run UMAP")
p1=proc.time()
sce <- runUMAP(sce, exprs_values=exprs,
    external_neighbors=TRUE, BPPARAM=bpParam,
    BNPARAM = knnParam,
    approx_pow=TRUE, n_threads=cores, pca=50, n_sgd_threads=cores)
p2=proc.time()
time_delta<-p2[3]-p1[1]
sprintf("run UMAP use %f seconds",time_delta)

print("### Run TSNE")
p1=proc.time()
sce <- runTSNE(sce, exprs_values=exprs,
    external_neighbors=TRUE, BPPARAM=bpParam,
    BNPARAM = knnParam,
    approx_pow=TRUE, n_threads=cores, pca=50, n_sgd_threads=cores)

p2=proc.time()
time_delta<-p2[3]-p1[1]
sprintf("run TSNE use %f seconds",time_delta)
# BNPARAM = KmknnParam()
#sce <- runUMAP(sce, exprs_values = "exprs")
#sce <- runTSNE(sce, exprs_values = "exprs")
saveRDS(sce,file.path(args$outdir,"sce_parallel.rds"))
