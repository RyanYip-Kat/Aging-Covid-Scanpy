library(Seurat)
library(argparse)
library(stringr)
library(ggplot2)

parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--outdir",
                    type="character",
                    default="output",
                    help="the path to save result")

args <- parser$parse_args()

if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}

my_theme<-theme(axis.title.x = element_text(size=25),
                  axis.text.x = element_text(size=18),
                  axis.text.y = element_text(size=18),
                  axis.title.y = element_text(size=25),
                  plot.title=element_text(size=25,face="bold"),
                  legend.text = element_text(size=15),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank())


seurat=readRDS("output/aging/NK/seurat.rds")
seurat=RenameIdents(seurat,"9"="NK1","15"="NK1","1"="NK2","4"="NK2","5"="NK2","13"="NK2",
		    "3"="NK3","6"="NK3","11"="NK3","12"="NK3","0"="NK3",
		    "7"="NK3","10"="NK3","8"="NK3","14"="NK3","2"="NK3")

seurat$celltype=Idents(seurat)
pdf(file.path(args$outdir,"celltype.pdf"),width=12,height=12)
DimPlot(seurat,reduction="tsne",group.by="celltype",pt.size=1.5)+my_theme
dev.off()

mat=as.matrix(table(seurat$sample_id,seurat$celltype))
write.table(mat,file.path(args$outdir,"ident_celltype_number.csv"),sep=",",quote=F)

mat=as.matrix(table(seurat$condition,seurat$celltype))
write.table(mat,file.path(args$outdir,"condition_celltype_number.csv"),sep=",",quote=F)


