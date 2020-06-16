library(CATALYST)
library(cowplot)
library(flowCore)
library(diffcyt)
library(scater)
library(SingleCellExperiment)
# R version 3.6.3 (2020-02-29)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 18.04.3 LTS

# Matrix products: default
# BLAS:   /usr/lib/x86_64-linux-gnu/openblas/libblas.so.3
# LAPACK: /usr/lib/x86_64-linux-gnu/libopenblasp-r0.2.20.so

# locale:
#  [1] LC_CTYPE=zh_CN.UTF-8       LC_NUMERIC=C
#  [3] LC_TIME=zh_CN.UTF-8        LC_COLLATE=en_US.UTF-8
#  [5] LC_MONETARY=zh_CN.UTF-8    LC_MESSAGES=en_US.UTF-8
#  [7] LC_PAPER=zh_CN.UTF-8       LC_NAME=C
#  [9] LC_ADDRESS=C               LC_TELEPHONE=C
# [11] LC_MEASUREMENT=zh_CN.UTF-8 LC_IDENTIFICATION=C

# attached base packages:
# [1] parallel  stats4    stats     graphics  grDevices utils     datasets
# [8] methods   base

# other attached packages:
#  [1] scater_1.14.6               ggplot2_3.3.0
#  [3] SingleCellExperiment_1.8.0  SummarizedExperiment_1.16.1
#  [5] DelayedArray_0.12.2         BiocParallel_1.20.1
#  [7] matrixStats_0.56.0          Biobase_2.46.0
#  [9] GenomicRanges_1.38.0        GenomeInfoDb_1.22.0
# [11] IRanges_2.20.2              S4Vectors_0.24.3
# [13] BiocGenerics_0.32.0         flowCore_1.52.1
# [15] cowplot_1.0.0               CATALYST_1.10.1

# loaded via a namespace (and not attached):
#   [1] readxl_1.3.1                circlize_0.4.8
#   [3] drc_3.0-1                   plyr_1.8.6
#   [5] igraph_1.2.5                lazyeval_0.2.2
#   [7] ConsensusClusterPlus_1.50.0 shinydashboard_0.7.1
#   [9] splines_3.6.3               fda_2.4.8.1
#  [11] TH.data_1.0-10              digest_0.6.25
#  [13] htmltools_0.4.0             viridis_0.5.1
#  [15] magrittr_1.5                CytoML_1.12.1
#  [17] cluster_2.1.0               ks_1.11.7
#  [19] openxlsx_4.1.4              limma_3.42.2
#  [21] ComplexHeatmap_2.2.0        RcppParallel_5.0.0
#  [23] R.utils_2.9.2               sandwich_2.5-1
#  [25] flowWorkspace_3.34.1        jpeg_0.1-8.1
#  [27] colorspace_1.4-1            rrcov_1.5-2
#  [29] ggrepel_0.8.2               haven_2.2.0
#  [31] dplyr_0.8.5                 crayon_1.3.4
#  [33] RCurl_1.98-1.1              jsonlite_1.6.1
#  [35] hexbin_1.28.1               graph_1.64.0
#  [37] survival_3.1-11             zoo_1.8-7
#  [39] glue_1.3.2                  flowClust_3.24.0
#  [41] gtable_0.3.0                nnls_1.4
#  [43] zlibbioc_1.32.0             XVector_0.26.0
#  [45] GetoptLong_0.1.8            ggcyto_1.14.1
#  [47] BiocSingular_1.2.2          car_3.0-7
#  [49] IDPmisc_1.1.20              Rgraphviz_2.30.0
#  [51] shape_1.4.4                 DEoptimR_1.0-8
#  [53] abind_1.4-5                 scales_1.1.0
#  [55] mvtnorm_1.1-0               Rcpp_1.0.4
#  [57] plotrix_3.7-7               xtable_1.8-4
#  [59] viridisLite_0.3.0           clue_0.3-57
#  [61] rsvd_1.0.3                  openCyto_1.24.0
#  [63] foreign_0.8-76              mclust_5.4.5
#  [65] FlowSOM_1.18.0              tsne_0.1-3
#  [67] DT_0.13                     httr_1.4.1
#  [69] htmlwidgets_1.5.1           RColorBrewer_1.1-2
#  [71] pkgconfig_2.0.3             XML_3.99-0.3
#  [73] R.methodsS3_1.8.0           flowViz_1.50.0
#  [75] reshape2_1.4.3              flowStats_3.44.0
#  [77] tidyselect_1.0.0            rlang_0.4.5
#  [79] later_1.0.0                 munsell_0.5.0
#  [81] cellranger_1.1.0            tools_3.6.3
#  [83] ggridges_0.5.2              shinyBS_0.61
#  [85] fastmap_1.0.1               stringr_1.4.0
#  [87] yaml_2.2.1                  zip_2.0.4
#  [89] robustbase_0.93-6           purrr_0.3.3
#  [91] RBGL_1.62.1                 mime_0.9
#  [93] R.oo_1.23.0                 compiler_3.6.3
#  [95] beeswarm_0.2.3              plotly_4.9.2
#  [97] curl_4.3                    png_0.1-7
#  [99] tibble_2.1.3                pcaPP_1.9-73
# [101] stringi_1.4.6               forcats_0.5.0
# [103] lattice_0.20-40             Matrix_1.2-18
# [105] shinyjs_1.1                 vctrs_0.2.4
# [107] pillar_1.4.3                lifecycle_0.2.0
# [109] GlobalOptions_0.1.1         BiocNeighbors_1.4.2
# [111] irlba_2.3.3                 data.table_1.12.8
# [113] bitops_1.0-6                corpcor_1.6.9
# [115] httpuv_1.5.2                R6_2.4.1
# [117] latticeExtra_0.6-29         promises_1.1.0
# [119] KernSmooth_2.23-16          gridExtra_2.3
# [121] rio_0.5.16                  vipor_0.4.5
# [123] codetools_0.2-16            MASS_7.3-51.5
# [125] gtools_3.8.1                assertthat_0.2.1
# [127] rjson_0.2.20                withr_2.1.2
# [129] mnormt_1.5-6                multcomp_1.4-12
# [131] GenomeInfoDbData_1.2.2      hms_0.5.3
# [133] ncdfFlow_2.32.0             grid_3.6.3
# [135] tidyr_1.0.2                 DelayedMatrixStats_1.8.0
# [137] carData_3.0-3               Rtsne_0.15
# [139] shiny_1.4.0.2               base64enc_0.1-3
# [141] ggbeeswarm_0.6.0            ellipse_0.4.1


safeBPParam <- function(nworkers) {
    if (.Platform$OS.type=="windows") {
        BiocParallel::SerialParam()
    } else {
        BiocParallel::MulticoreParam(nworkers)
    }
}
(cores <- parallel::detectCores())
(bpParam <- safeBPParam(cores))


files=file.path("output/cutoff/data",list.files("output/cutoff/data",pattern=".fcs"))
pbmc_fs <- read.flowSet(files, transformation = FALSE, truncate_max_range = FALSE)

pbmc_panel<-readRDS("output/cutoff/data/pbmc_panel.rds")
pbmc_md<-readRDS("output/cutoff/data/pbmc_md.rds")

(sce <- prepData(pbmc_fs, pbmc_panel, pbmc_md))

# options(verbose=TRUE)

# specify markers to use for clustering
lineage_markers <- rownames(sce)
sce <- cluster(sce, features = lineage_markers, xdim = 5, ydim = 5, maxK = 10, verbose = FALSE, seed = 1)

exprs <- "exprs"
# # KNN 的计算在 UMAP 和 TSNE 是可以复用的
# library(BiocNeighbors)
# x <- assay(sce, exprs_values = exprs)
# x_index <- buildKmknn(x)
# (knnParam <- KmknnParam(BNINDEX=x_index))

# from https://github.com/jlmelville/uwot
# For umap, it's better to provide a and b directly with a fixed precision rather than allowing them to be calculated via the spread and min_dist parameters. For default UMAP, use a = 1.8956, b = 0.8006.
# Using n_sgd_threads with more than 1 thread will not give reproducible results, but should not behave any worse than LargeVis in that regard, so for many visualization needs, this is also worth trying.
# You can (and should) adjust the number of threads via the n_threads, which controls the nearest neighbor and smooth knn calibration, and the n_sgd_threads parameter, which controls the number of threads used during optimization. For the n_threads, the default is half of whatever RcppParallel thinks should be the default. For n_sgd_threads the default is 0, which ensures reproducibility of results with a fixed seed.
# Rprof("Rprof-runUMAP.out")
system.time(
    sce <- runUMAP(sce, exprs_values=exprs, 
    external_neighbors=TRUE, BPPARAM=bpParam, 
    approx_pow=TRUE, init="spca", n_threads=cores, pca=50, n_sgd_threads=cores)
)
# Rprof(NULL)