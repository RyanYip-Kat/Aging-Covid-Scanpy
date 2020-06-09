library(argparse)
library(harmony)
library(Seurat)
library(stringr)
library(sva)

parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--outdir",
                    type="character",
                    default="output",
                    help="the path to save result")

parser$add_argument("--seurat",
                    type="character",
                    default="")

args <- parser$parse_args()

if(!dir.exists(args$outdir)){o
        dir.create(args$outdir,recursive=TRUE)
}

seurat<-readRDS(args$seurat)
Names<-colnames(seurat)
#status=unlist(lapply(Names,function(x){return(str_split(x,"-")[[1]][1])}))
sample=unlist(lapply(Names,function(name){return(str_split(name,"-")[[1]][2])}))

metadata<-data.frame("cell_id"=Names,"sample"=sample) #"status"=status
print("### RunHarmony")
if("pca"%in%names(seurat@reductions)){
	mat=Embeddings(seurat,"pca")[,1:20]
        harmony_embeddings<-HarmonyMatrix(mat, meta_data=metadata, vars_use=c("sample"),do_pca=FALSE) #v3
}else{
	mat=GetAssayData(seurat,"counts")
	harmony_embeddings<-HarmonyMatrix(mat, meta_data=metadata, vars_use=c("sample"),do_pca=TRUE)
}

rownames(harmony_embeddings)<-Names
colnames(harmony_embeddings)=paste("Harmony",1:ncol(harmony_embeddings),sep="_")
saveRDS(harmony_embeddings,file.path(args$outdir,"harmony.rds"))

print("### Add harmony")
mat<-as.matrix(harmony_embeddings)
colnames(mat)<-paste("Harmony_",1:ncol(mat),sep = "")
seurat[["harmony"]]<-CreateDimReducObject(embeddings =mat,
                                    key = "Harmony_",
                                    assay = DefaultAssay(seurat))
print("### Run TSNE")
seurat <- RunTSNE(object = seurat, reduction = "harmony", dims = 1:20)
print("Run FindNeighbors")
seurat <- FindNeighbors(object = seurat,reduction = "harmony", dims = 1:20)

print("Run Find Clusters")
seurat <- FindClusters(object = seurat,resolution=0.8, verbose =TRUE)


saveRDS(seurat,file.path(args$outdir,"seurat_harmony.rds"))

