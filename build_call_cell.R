library(argparse)
library(stringr)
library(readr)

parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--outdir",
                    type="character",
                    default="output",
                    help="the path to save result")

parser$add_argument("--path",
                    type="character",
                    default="")

args <- parser$parse_args()

if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}

peak=file.path(args$path,"peaks.bed")
barcode=file.path(args$path,"barcodes.tsv")

peak=read.table(peak)
peak$V1=as.character(peak$V1)
barcode=read.table(barcode)
write.table(barcode,file.path(args$outdir,"barcodes.txt"),sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)

peak_list=list()
for(i in 1:nrow(peak)){
	x<-peak[i,]
	p<-paste(x,collapse="-")
	peak_list[[i]]=p
}

df=as.data.frame(do.call(rbind,peak_list))
write.table(df,file.path(args$outdir,"features.txt"),sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)

