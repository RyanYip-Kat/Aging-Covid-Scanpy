library(argparse)
library(stringr)
library(chromVAR)
library(BiocParallel)
library(RColorBrewer)
library(viridis)
library(pheatmap)
library(data.table)
GENOME_NAME = 'hg38'
#############################
# Do DA with one cluster vs the rest clusters
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
#############################
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--down_dir",
                    type="character",
                    default="")

parser$add_argument("--cluster",
                    type="character",
                    default=NULL)

args <- parser$parse_args()


down_dir=args$down_dir
seurat_obj<-readRDS(file.path(down_dir, 'seurat_obj.rds'))
metaData = seurat_obj@meta.data
head(metaData)

chromVar.obj = readRDS(file.path(down_dir, 'chromVar_obj.rds'))

if(is.null(args$cluster)){
	motif_enriched="differential_TF_motif_enriched_in_clusters.txt"
}else{
	motif_enriched=paste0(args$cluster,"_","differential_TF_motif_enriched_in_clusters.txt")
}

diff_tf_enrich_file = file.path(down_dir, motif_enriched)

if(!file.exists(diff_tf_enrich_file)){
  dev = deviations(chromVar.obj)
  da.res = do_DA_motif(dev,
                       clusters = data.table('barcode' = rownames(metaData),
                                             'cluster' = metaData[[args$cluster]]),
                       topn = 10)
  rm(dev)
}else{
  da.res = fread( paste0(down_dir, '/differential_TF_motif_enriched_in_clusters.txt'))
}

## plot enriched TFs in heatmap
sele.tfs = da.res$feature
zscores = deviationScores(chromVar.obj)
sele.zscores = zscores[sele.tfs, ]


## change tf name to be more readable
if(grepl(GENOME_NAME, pattern = 'hg', ignore.case = T)){
  rnames = rownames(sele.zscores)
  nnames = sapply(rnames, function(x) unlist(strsplit(x, '_'))[3])
  nnames1 = sapply(rnames, function(x) unlist(strsplit(x, '_'))[1])
  rownames(sele.zscores) = ifelse(grepl(nnames, pattern = 'LINE'), nnames1, nnames)
}else{
  rnames = rownames(sele.zscores)
  nnames = sapply(rnames, function(x) unlist(strsplit(x, '_'))[3])
  rownames(sele.zscores) = nnames
  sele.zscores = sele.zscores[!grepl(nnames, pattern = '^LINE'), ]
}


#metaData$active_clusters = as.character(metaData$active_clusters)
metaData$active_clusters = as.character(metaData[[args$cluster]])
metaData$active_clusters<-factor(metaData$active_clusters,levels=c("NKT","BC","MC","undefined"))
metaData = data.table(metaData, keep.rownames = T)
setkey(metaData, active_clusters)
head(metaData)

rr = metaData$rn[metaData$rn %in% colnames(sele.zscores)]
sele.zscores = sele.zscores[, rr]


sele.zscores = sele.zscores[!duplicated(sele.zscores), ]

ann_col = data.frame('cluster' = metaData$active_clusters)
rownames(ann_col) = metaData$rn

up_cut = quantile(sele.zscores, 0.95, na.rm = T)
low_cut = quantile(sele.zscores, 0.05, na.rm = T)
sele.zscores[is.na(sele.zscores)] = 0
low_cut = min(0, low_cut)
sele.zscores[sele.zscores > up_cut] = up_cut
sele.zscores[sele.zscores < low_cut] = low_cut

#cluster = brewer.pal(n=length(unique(metaData$active_clusters)), name = 'Paired')
#names(cluster) = sort(unique(metaData$active_clusters))
#ann_colors = list('cluster' = cluster)

# resample to reduce memory used
set.seed(2019)
rids = sort(sample(1:ncol(sele.zscores), floor(ncol(sele.zscores)/6)))
ann_col0 = data.frame(ann_col[rids, ])
rownames(ann_col0) = colnames(sele.zscores)[rids]
mtx0 = sele.zscores[, rids]
names(ann_col0) = 'cluster'

print("### heatmap")
pdf(file.path(down_dir,paste0(args$cluster,"_heatmap_motif_enrich.pdf")),width=12,height=12)
#pheatmap::pheatmap(mtx0, cluster_cols = F,
#                         cluster_rows = T, show_colnames = F, fontsize = 13,
#                         annotation_col = ann_col0, color = viridis(100),
#                         annotation_colors = ann_colors, fontsize_row = 9)

pheatmap::pheatmap(mtx0, cluster_cols = F,
		   cluster_rows = T, show_colnames = F, fontsize = 13,
		   annotation_col = ann_col0, color = viridis(100),fontsize_row = 10)
dev.off()
