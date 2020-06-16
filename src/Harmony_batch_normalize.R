library(harmony)
library(argparse)
library(stringr)
library(sva)
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--count",
                    type="character",
                    default="",
                    help="the path of fcs files")


parser$add_argument("--outdir",
                      type="character",
                      default="output",
                      help="save path")

args<-parser$parse_args()
if(!dir.exists(args$outdir)){
            dir.create(args$outdir,recursive=TRUE)
}

counts<-readRDS(args$count)                      
Names<-rownames(counts)

status=unlist(lapply(Names,function(x){return(str_split(x,"-")[[1]][1])}))
sample=unlist(lapply(Names,function(name){return(str_split(name,"_")[[1]][2])}))
metadata<-data.frame("cell_id"=Names,"sample"=sample,"status"=status)
harmony_embeddings<-HarmonyMatrix(counts, metadata, c("status"))
rownames(harmony_embeddings)<-Names
colnames(harmony_embeddings)=paste("Harmony",1:ncol(harmony_embeddings),sep="_")
saveRDS(harmony_embeddings,file.path(args$outdir,"harmony.rds"))

