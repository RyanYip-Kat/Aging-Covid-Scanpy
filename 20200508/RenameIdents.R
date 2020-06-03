library(Seurat)
seurat=readRDS("output/aging/20000/seurat.rds")
seurat=RenameIdents(seurat,"0"="TC","6"="TC","8"="TC","17"="TC",
		    "24"="TC","10"="TC","13"="TC","12"="TC","9"="TC","11"="TC","4"="TC",
		    "16"="DC","22"="DC","18"="DC","20"="DC",
		    "5"="BC","14"="BC","25"="BC","26"="BC",
		    "1"="Mono","3"="Mono","21"="Mono",
		    "2"="NK","7"="NK","15"="NK","19"="NK",
		    "23"="Unknow")

print(table(Idents(seurat)))

pdf("celltype.pdf",width=12,height=12)
DimPlot(seurat,reduction="tsne")
dev.off()
