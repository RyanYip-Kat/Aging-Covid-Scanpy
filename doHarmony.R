library(harmony)
library(Seurat)
library(stringr)
library(sva)

runHarmony<-function(seurat,npcs=20,sample_id=NULL){
	Names<-colnames(seurat)
	if(is.null(sample_id)){
		sample_id=unlist(lapply(Names,function(name){return(str_split(name,"-")[[1]][2])}))
	}
	metadata<-data.frame("cell_id"=Names,"sample"=sample_id) #"status"=status
	print("### RunHarmony")
	stopifnot("pca"%in%names(seurat@reductions)
	mat=Embeddings(seurat,"pca")[,1:npcs]
	harmony_embeddings<-HarmonyMatrix(mat, meta_data=metadata, vars_use=c("sample"),do_pca=FALSE) #v3
	rownames(harmony_embeddings)<-Names
	colnames(harmony_embeddings)=paste("Harmony",1:ncol(harmony_embeddings),sep="_")

	print("### Add harmony")
        mat<-as.matrix(harmony_embeddings)
        colnames(mat)<-paste("Harmony_",1:ncol(mat),sep = "")
        seurat[["harmony"]]<-CreateDimReducObject(embeddings =mat,
                                    key = "Harmony_",
                                    assay = DefaultAssay(seurat))
        print("### Run TSNE")
        seurat <- RunTSNE(object = seurat, reduction = "harmony", dims = 1:20)
        print("### Run UMAP")
	seurat <- RunUMAP(object = seurat, reduction = "harmony", dims = 1:20)
        print("Run FindNeighbors")
        seurat <- FindNeighbors(object = seurat,reduction = "harmony", dims = 1:20)

        print("Run Find Clusters")
        seurat <- FindClusters(object = seurat,resolution=0.8, verbose =TRUE)
	return(seurat)
}


