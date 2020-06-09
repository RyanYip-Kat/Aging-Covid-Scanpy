library(parallel)
library(chromVAR)
library(Seurat)
library(argparse)
library(data.table)
## plot heatmap ####
library(BiocParallel)
register(SerialParam())

# Do DA/DE with one cluster vs the rest clusters
# clusters are the data frame with <barcode> <cluster>
do_DA_motif <- function(mtx_score, clusters, test = 'wilcox', 
                  only.pos = T, fdr = 0.05, topn = 10){
  clusters$cluster = as.character(clusters$cluster)
  cls = unique(clusters$cluster)
  res = NULL
  features = rownames(mtx_score)
  for(cluster0 in cls){
    bc0 = clusters[cluster == cluster0]$barcode
    mtx1 = mtx_score[, colnames(mtx_score) %in% bc0]
    mtx2 = mtx_score[, !colnames(mtx_score) %in% bc0]
    mu1 = sapply(1:length(features), function(x) mean(mtx1[x, ]))
    mu2 = sapply(1:length(features), function(x) mean(mtx2[x, ]))
    
    pvs = rep(0.5, length(features))
    
    for(x in 1:length(features)){
      a1 = mtx1[x, ]
      a2 = mtx2[x, ]
      if(length(which(!is.na(a1))) < 2 || length(which(!is.na(a2))) < 2) next
      pvs[x] = wilcox.test(a1, a2, alternative = 'greater')$p.value
    }
    
    pvs.adj = p.adjust(pvs, method = 'fdr')
    res0 = data.table('feature' = features, 'cluster' = cluster0,
                      'mean1' = mu1, 'mean2' = mu2,
                      'pv' = pvs, 'pv_adjust' = pvs.adj)
    
    
    res0 = res0[order(pv_adjust), ]
    res0 = res0[pv_adjust <= fdr]
    
    if(nrow(res0) > topn) res0 = res0[1:topn, ]
    res = rbind(res, res0)
  }
  return(res)
}
###########################
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--output_dir",
                    type="character",
                    default="")

parser$add_argument("--active_clusters",
                    type="character",
                    default="seurat_clusters")

args <- parser$parse_args()
output_dir=args$output_dir

seurat_file <- paste0(output_dir, '/seurat_obj.rds')
ss = readRDS(seurat_file)
chromVar.obj = readRDS(paste0(output_dir, '/chromVar_obj.rds'))

metaData = ss@meta.data

dev = deviations(chromVar.obj)
da.res = do_DA_motif(dev, 
                 clusters = data.table('barcode' = rownames(metaData),
                                       'cluster' = metaData[[args$active_clusters]]),topn = 10)
write.csv(da.res, file = file.path(output_dir, paste0(args$active_clusters,'_differential_TF_motif_enriched_in_clusters.txt')), 
              quote = F, row.names = F )
    
    
