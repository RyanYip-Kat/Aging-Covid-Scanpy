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


seurat=readRDS("../20200529/output/aging/TC/seurat.rds")
Idents(seurat)=seurat$seurat_clusters
seurat=RenameIdents(seurat,"10"="CD4+CD8+","26"="CD4+CD8+","4"="CD4-CD8-",
		    "6"="CD4 Naive","7"="CD4 Naive","22"="CD4 Naive","24"="CD4 Naive",
		    "0"="CD4 Tcm","2"="CD4 Tem","5"="CD4 Tem","8"="CD4 Tem","14"="CD4 Tem",
		    "16"="CD4 Tem","19"="CD4 Treg","1"="CD8 Naive","3"="CD8 Naive","20"="CD8 Naive",
		    "21"="CD8 Naive","27"="CD8 Naive","9"="CD8 Tem","11"="CD8 Tem","12"="CD8 Tem","13"="CD8 Tem",
		    "17"="CD8 Tem","25"="CD8 Tem","15"="CD8 CTL","18"="CD8 CTL","23"="CD8 CTL")

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

