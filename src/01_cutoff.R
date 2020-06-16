library(CATALYST)
library(stringr)
library(argparse)
library(flowCore)
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--path",
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

DATA<-read.csv("20200330/result/panel.csv",stringsAsFactors=F)
name=DATA$name
fcs_path<-file.path(args$path,list.files(args$path,pattern=".fcs"))
for(file in fcs_path){
	ff<-flowCore::read.FCS(file)
	data<-ff@parameters@data
	print(paste0("Number of protein : ",length(name)))
	ff@parameters@data<-data[data$name%in%name,]
	ff@exprs<-ff@exprs[,name]

	write.FCS(ff,file.path(args$outdir,basename(file)))
}
