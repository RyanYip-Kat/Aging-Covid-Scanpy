library(argparse)
library(stringr)
library(pheatmap)

parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--outdir",
                    type="character",
                    default="output",
                    help="the path to save result")

parser$add_argument("--down_dir",
                    type="character",
                    default="")

args <- parser$parse_args()

if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}

down_dir=args$down_dir

group1_fp = '0'
group2_fp = 'rest'
footprint_stats.file = paste0(down_dir, '/differential_TF_footprint_', 
                              group1_fp, '_vs_', group2_fp, '.txt')
if(file.exists(footprint_stats.file)){
  
  footprint_out = fread(footprint_stats.file)
  if(length(unique(footprint_out$motif)) > 100){
  footprint_out[, 'N' := .N, by = higher_in_cluster]
  cls = unique(footprint_out[N > 10]$higher_in_cluster)
  if(length(cls) >= 1){
    res0 = NULL
    for(cl0 in cls){
      tmp = footprint_out[higher_in_cluster == cl0]
      tmp = tmp[order(P_values)][1:10, ]
      res0 = rbind(res0, tmp)
    }
    footprint_out = rbind(footprint_out[N < 10], res0)
  }
}

  mm = reshape2::acast(motif ~ higher_in_cluster, data = footprint_out, 
                       value.var = "P_values")
  mm = -log10(mm)
  mm[is.na(mm)] = 0
  cn = colnames(mm)
  cn.new = sapply(cn, function(x) gsub('_higher', '', x))
  colnames(mm) = cn.new
  mm[mm > 3] = 3

  pdf(file.path(args$outdir,paste0(group1_fp,'_vs_', group2_fp,".pdf")),width=10,height=16)
  pheatmap(mm, cluster_cols = F, fontsize = 13, fontsize_row = 9,
                           color = viridis::viridis(100))
  
  dev.off()

}
