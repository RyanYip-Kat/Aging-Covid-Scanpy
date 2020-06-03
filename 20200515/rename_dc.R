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


seurat=readRDS("../20200514/new.Cluster/aging/DC/seurat.rds")
seurat=RenameIdents(seurat,"3"="pDC","4"="pDC","5"="pDC","7"="pDC","9"="pDC",
		    "10"="pDC","11"="pDC","12"="pDC","16"="pDC","20"="pDC",
		    "0"="cDC","1"="cDC","2"="cDC","6"="cDC","8"="cDC","13"="cDC",
		    "14"="cDC","15"="cDC","17"="cDC","18"="cDC","19"="cDC")

seurat$celltype=Idents(seurat)
print(table(seurat$celltype))
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
