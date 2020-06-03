library(Seurat)
library(argparse)
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--outdir",
                    type="character",
                    default="output",
                    help="the path to save result")

args <- parser$parse_args()

if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}

seurat=readRDS("new.Cluster/aging/Total/seurat.rds")
Idents(seurat)=seurat$new.Clusters
seurat=RenameIdents(seurat,"1"="BC","14"="BC","29"="BC","34"="BC","31"="BC",
		    "18"="DC","33"="DC","22"="DC","26"="DC",
		    "15"="NK","32"="NK","27"="NK","5"="NK","4"="NK","23"="NK",
		    "6"="Mono","7"="Mono","9"="Mono","10"="Mono","25"="Mono","30"="Mono",
		    "24"="TC","2"="TC","19"="TC","16"="TC","17"="TC","20"="TC",
		    "12"="TC","11"="TC","13"="TC","0"="TC","3"="TC","8"="TC",
		    "21"="TC","35"="TC","28"="Basophil","TC"="TC","NK"="NK","Mono"="Mono")

print(table(Idents(seurat)))
seurat$celltype=Idents(seurat)
pdf(file.path(args$outdir,"celltype.pdf"),width=12,height=12)
DimPlot(seurat,reduction="tsne")
dev.off()

mat=table(seurat$sample_id,seurat$celltype)
write.table(as.matrix(mat),file.path(args$outdir,"ident_cluster_numbers.csv"),sep=",",quote=F)

mat=table(seurat$celltype)
write.table(as.matrix(mat),file.path(args$outdir,"cluster_numbers.csv"),sep=",",quote=F)

saveRDS(seurat,file.path(args$outdir,"seurat.rds"))
