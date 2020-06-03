library(Seurat)
library(pheatmap)

#seurat=readRDS("../20200602/output/COVID-CYTOF/Total/seurat.rds")
seurat=readRDS("../20200529/output/aging/Total/seurat.rds")

path="aging-clusters"
files=list.files(path,".csv",recursive=T,full.names=TRUE)
DATA=lapply(files,function(file){return(read.csv(file,stringsAsFactors=F))})
DATA=do.call(rbind,DATA)
rownames(DATA)=DATA$barcode
DATA=DATA[,"celltype",drop=F]
seurat<-AddMetaData(seurat,metadata=DATA,col.name="celltype")
Idents(seurat)=seurat$celltype

markers=c("CD3","CD4","CD8A","CD45RO","CD45RA","CCR7","CCR6","CD127","CD27","CD161","CXCR3","CD56","CD57","CD19","CD20", "IGD","CXCR5","CD38","CD11C","CD123","HLA-DR","CD25","CCR4","CD14","CD16","CD294")
cts=AverageExpression(seurat,features=markers)[[DefaultAssay(seurat)]]

#cts=as.data.frame(cts)
#cts=as.matrix(cts)
pdf("output/aging/Total/subset_block_heatmap1.pdf",width=10,height=16)
p=pheatmap(cts,
            border_color=NA,
            show_rownames=T,
            show_colnames=T,
            scale="column",
            cluster_rows=F,
            cluster_cols=F,fontsize=10)

print(p)
dev.off()


pdf("output/aging/Total/subset_block_heatmap2.pdf",width=16,height=10)
p=pheatmap(t(cts),
            border_color=NA,
            show_rownames=T,
            show_colnames=T,
            scale="column",
            cluster_rows=F,
            cluster_cols=F,fontsize=10)

print(p)
dev.off()
