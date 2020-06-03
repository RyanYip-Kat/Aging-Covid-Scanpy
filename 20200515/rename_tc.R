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


seurat=readRDS("../20200514/new.Cluster/aging/TC/seurat.rds")
seurat=RenameIdents(seurat,"9"="CD4+CD8+","25"="CD4+CD8+",
		    "3"="CD4-CD8-",
		    "5"="CD4 Naive","8"="CD4 Naive",
		    "1"="CD4 Tcm","20"="CD4 Tcm",
		    "2"="CD4 Tem","4"="CD4 Tem","6"="CD4 Tem","15"="CD4 Tem",
		    "22"="CD4 Tem","23"="CD4 Tem","18"="CD4 Treg","17"="CD4 CTL",
		    "0"="CD8 Naive","13"="CD8 Naive","19"="CD8 Naive","24"="CD8 Naive",
		    "7"="CD8 Tem","10"="CD8 Tem","11"="CD8 Tem","16"="CD8 Tem","26"="CD8 Tem",
		    "12"="CD8 CTL","14"="CD8 CTL","21"="CD8 Tem")

seurat$celltype=Idents(seurat)
#seurat$celltype=factor(seurat$celltype,levels=c("CD4 Naive","CD4 Tcm","CD4 Tem","CD4 Treg","CD8 Naive","CD8 Tem","CD8 CTL","NKT","CD4+CD8+","CD4-CD8-"))
pdf(file.path(args$outdir,"celltype.pdf"),width=12,height=12)
DimPlot(seurat,reduction="tsne",group.by="celltype")+my_theme
dev.off()

mat=as.matrix(table(seurat$sample_id,seurat$celltype))
write.table(mat,file.path(args$outdir,"ident_celltype_number.csv"),sep=",",quote=F)

mat=as.matrix(table(seurat$condition,seurat$celltype))
write.table(mat,file.path(args$outdir,"condition_celltype_number.csv"),sep=",",quote=F)

mat=as.matrix(table(seurat$celltype))
write.table(mat,file.path(args$outdir,"celltype_number.csv"),sep=",",quote=F)

saveRDS(seurat,file.path(args$outdir,"seurat.rds"))

