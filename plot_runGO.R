library(argparse)
library(stringr)
library(ggplot2)
library(xlsx)

library(data.table)
library(magrittr)

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
group1 = '0'
group2 = '1'
go_file = paste0(down_dir, 'enrichedGO_differential_accessible_features_', group1, '_vs_', group2, '.xlsx')

# get enriched terms for cluster0 for example and show top 20 terms
go_res = xlsx::read.xlsx(go_file, sheetName = 'cluster0')

go_res = data.table(go_res)
go_res[, 'score' := -log10(p.adjust)]
go_res = go_res[order(-score), ]
ngo = min(20, nrow(go_res))
go_res = go_res[1:ngo, ]
go_res = go_res[order(score), ]
go_res$Description = factor(go_res$Description, levels = go_res$Description)
        
p_go <- ggplot(go_res, aes(y = score, x = Description, fill = Count)) +
  geom_bar(width = 0.7, stat = 'identity') +
  ggtitle("Enriched terms: cluster_0") + theme_classic() + 
  theme(legend.position = 'bottom', legend.direction = "horizontal") + 
  coord_flip()  +  scale_fill_continuous(name = "#genes", type = "viridis") +
  xlab('') + ylab('-log10(p.adjust)')
         
pdf(file.path(args$outdir,paste0("GO_",group1,"_vs_",group2,".pdf")),width=12,height=16)
print(p_go)
dev.off()
