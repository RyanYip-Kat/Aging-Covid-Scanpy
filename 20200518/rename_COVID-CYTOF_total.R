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


seurat=readRDS("output/COVID-CYTOF/20000/seurat.rds")
seurat=RenameIdents(seurat,"1"="TC","5"="TC","7"="TC",
		    "10"="TC","12"="TC","15"="TC","16"="TC","19"="TC","21"="TC",
		    "24"="TC","2"="NK","3"="NK","8"="NK","11"="NK",
		    "20"="NK","23"="NK","25"="NK","28"="NK","4"="BC","22"="BC","27"="BC",
		    "0"="Mono","6"="Mono","9"="Mono","17"="Mono","26"="Mono",
		    "14"="DC","18"="DC","13"="Basophil")


seurat$celltype=Idents(seurat)
pdf(file.path(args$outdir,"celltype.pdf"),width=12,height=12)
DimPlot(seurat,reduction="tsne",group.by="celltype")+my_theme
dev.off()

mat=as.matrix(table(seurat$sample_id,seurat$celltype))
write.table(mat,file.path(args$outdir,"ident_celltype_number.csv"),sep=",",quote=F)

mat=as.matrix(table(seurat$sample_id,seurat$seurat_clusters))
write.table(mat,file.path(args$outdir,"ident_clusters_number.csv"),sep=",",quote=F)


mat=as.matrix(table(seurat$condition,seurat$celltype))
write.table(mat,file.path(args$outdir,"condition_celltype_number.csv"),sep=",",quote=F)

mat=as.matrix(table(seurat$celltype))
write.table(mat,file.path(args$outdir,"celltype_number.csv"),sep=",",quote=F)

saveRDS(seurat,file.path(args$outdir,"seurat.rds"))
