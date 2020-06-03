library(Seurat)
library(argparse)
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


path="clusters"
files=list.files(path,".csv",recursive=T,full.names=TRUE)
DATA=lapply(files,function(file){return(read.csv(file,stringsAsFactors=F))})
DATA=do.call(rbind,DATA)
rownames(DATA)=DATA$barcode
DATA=DATA[,"celltype",drop=F]
seurat=readRDS("output/COVID-CYTOF/Total/seurat.rds")
seurat<-AddMetaData(seurat,metadata=DATA,col.name="celltype")
seurat$celltype=factor(seurat$celltype,levels=c("TC","NK","BC","Mono","DC"))

cols<-c("#1F77B4","#2DB8B9","#79BE5B","#E9726C","#F5A623")  # TC,NK,BC,Mono,DC

pdf(file.path(args$outdir,"celltype.pdf"),width=12,height=12)
DimPlot(seurat,reduction="tsne",group.by="celltype",cols=cols)+my_theme
dev.off()

mat=as.matrix(table(seurat$sample_id,seurat$celltype))
write.table(mat,file.path(args$outdir,"ident_celltype_number.csv"),sep=",",quote=F)

mat=as.matrix(table(seurat$condition,seurat$celltype))
write.table(mat,file.path(args$outdir,"condition_celltype_number.csv"),sep=",",quote=F)

mat=as.matrix(table(seurat$celltype))
write.table(mat,file.path(args$outdir,"celltype_number.csv"),sep=",",quote=F)


status<-with(seurat@meta.data,ifelse(condition%in%c("YH","AH"),"HC","CO"))
seurat$status<-factor(status,levels=c("HC","CO"))
pdf(file.path(args$outdir,"hc_co_celltype.pdf"),width=12,height=12)
DimPlot(seurat,reduction="tsne",group.by="celltype",split.by="status",cols=cols)+my_theme
dev.off()

Idents(seurat)=seurat$condition
seurat<-subset(seurat,idents=c("YCO","ACO"))
seurat$condition<-factor(seurat$condition,levels=c("YCO","ACO"))
pdf(file.path(args$outdir,"yco_aco_celltype.pdf"),width=12,height=12)
DimPlot(seurat,reduction="tsne",group.by="celltype",split.by="condition",cols=cols)+my_theme
dev.off()
