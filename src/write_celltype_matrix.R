library(Seurat)
library(stringr)
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

seurat=readRDS("../output/cohort/cohort-2/25000-V2/seurat.rds")
celltypes<-as.character(unique(Idents(seurat)))
for(celltype in celltypes){
	print(paste0("Loading ",celltype))
	x<-subset(seurat,idents=celltype)
	cts<-as.matrix(GetAssayData(x,slot="data"))
	#name=paste("Cell Line:",colnames(cts),sep=" ")
	#type=paste("CellType :",celltype,sep=" ")
	#status=unlist(lapply(colnames(cts),function(name){return(str_split(name,"-")[[1]][1])}))
	#status=paste("Status:",status,sep=" ")
	#DATA<-rbind(name,type,status,cts)
	#rownames(DATA)<-c("","","",paste("Gene:",rownames(cts),sep=" "))
        DATA<-as.matrix(cts)
	DATA<-t(DATA)
	filename<-file.path(args$outdir,paste0(celltype,".txt"))
	write.table(DATA,filename,sep="\t",quote=FALSE)
}



