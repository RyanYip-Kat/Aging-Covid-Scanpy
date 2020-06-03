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


seurat=readRDS("output/COVID-CYTOF/20000/NK/seurat.rds")
seurat=RenameIdents(seurat,"19"="NK1","30"="NK1",
		    "0"="NK2","4"="NK2","6"="NK2","7"="NK2","8"="NK2",
		    "10"="NK2","11"="NK2","14"="NK2","16"="NK2",
		    "17"="NK2","20"="NK2","22"="NK2","23"="NK2","24"="NK2","27"="NK2",
		    "29"="NK2","33"="NK2","1"="NK3","2"="NK3","3"="NK3","5"="NK3",
		    "9"="NK3","12"="NK3","13"="NK3","15"="NK3","18"="NK3","21"="NK3",
		    "25"="NK3","26"="NK3","28"="NK3","31"="NK3","32"="NK3")

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
