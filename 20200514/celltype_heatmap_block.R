library(Seurat)
library(pheatmap)

markers=c("CD3","CD4","CD8A","CD45RO","CD45RA","CCR7","CCR6","CD127","CD27","CD161","CXCR3","CD56","CD57","CD19","CD20", "IGD","CXCR5","CD38","CD11C","CD123","HLA-DR","CD25","CCR4","CD14","CD16","CD294")
seurat=readRDS("new.Cluster/aging/Total/seurat.rds")
cts=AverageExpression(seurat,features=markers)[[DefaultAssay(seurat)]]
#cts=as.data.frame(cts)
cts=cts[,c("TC","NK","BC","DC","Mono")]
#cts=as.matrix(cts)
pdf("new.Cluster/aging/Total/block_heatmap1.pdf",width=10,height=16)
p=pheatmap(cts,
            border_color=NA,
            show_rownames=T,
            show_colnames=T,
            scale="column",
            cluster_rows=F,
            cluster_cols=F,fontsize=10)

print(p)
dev.off()


pdf("new.Cluster/aging/Total/block_heatmap2.pdf",width=16,height=10)
p=pheatmap(t(cts),
            border_color=NA,
            show_rownames=T,
            show_colnames=T,
            scale="column",
            cluster_rows=F,
            cluster_cols=F,fontsize=10)

print(p)
dev.off()
