library(ggplot2)
library(argparse)
library(Seurat)
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--outdir",
                    type="character",
                    default="output",
                    help="the path to save result")


parser$add_argument("--seurat",
                    type="character",
                    default="")

parser$add_argument("--pt_size",
                    type="double",
                    default="1.0")

args <- parser$parse_args()
if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}

my_theme<-theme(axis.title.x = element_text(size=25),
                  axis.text.x = element_text(size=18),
                  axis.text.y = element_text(size=18),
                  axis.title.y = element_text(size=25),
                  plot.title=element_text(size=25,face="bold"),
                  legend.text = element_text(size=25),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank())

# cols = c("#191970","yellow","#FF3030")
seurat<-readRDS(args$seurat)
print(table(seurat$celltype))

pdf(file.path(args$outdir,"cluster.pdf"),width=12,height=12)
DimPlot(seurat,reduction="tsne",group.by="celltype",pt.size=args$pt_size)+my_theme
dev.off()

