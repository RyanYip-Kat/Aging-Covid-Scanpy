library(flowCore)
library(argparse)

parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--path",
                    type="character",
                    default="",
                    help="the path of fcs files")


args<-parser$parse_args()

path=args$path
files=file.path(path,list.files(path,".fcs"))

n<-0
for(file in files){
	x=read.FCS(file)
	n_row<-nrow(exprs(x))
	print(paste0("File : ",basename(file),"  number of cells :",n_row))
	n<-n+n_row
}

print(paste0("Total number :",n))
