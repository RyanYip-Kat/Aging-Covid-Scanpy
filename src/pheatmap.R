library(Seurat)
library(pheatmap)
color=colorRampPalette(c("navy","#40E0D0","#B0171F"))(50)
pheatmap.idents<-function(seurat,features){
        cts <- GetAssayData(seurat, slot = "data")
        features<-features[features%in%rownames(cts)]
        new_cluster <- sort(seurat$celltype)
        cts <- as.matrix(cts[features, names(new_cluster)])
        ac=data.frame(cluster=new_cluster)
        #color<-c("#40E0D0","#B0171F")
        p<-pheatmap(cts,color=color,show_colnames =F,show_rownames = T,
                    cluster_rows = F,
                    cluster_cols = F,
                    fontsize=10,
                    annotation_col=ac)
        return(p)
}

