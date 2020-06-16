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

fcs_path<-file.path(args$path,list.files(args$path,pattern="\\d+.fcs$"))
print(paste0("There are : ", length(fcs_path)," fcs files"))
ff<-concatFCS(fcs_path)
write.FCS(ff,file.path(args$outdir,"clean.fcs"))
