library(CATALYST)
library(stringr)
library(argparse)
library(flowCore)

parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--outdir",
                    type="character",
                    default="output",
                    help="the path to save result")


parser$add_argument("--fcs",
                    type="character",
                    default="")

args <- parser$parse_args()
if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}

ff<-flowCore::read.FCS(args$fcs,transformation = FALSE, truncate_max_range = FALSE)
DATA <- flowCore::exprs(ff)
print(paste0("data size : (",nrow(DATA),",",ncol(DATA),")"))

DATA<-as.data.frame(DATA)

RuleForCellType<-function(DATA){
	DATA<-subset(DATA,Ir193Di>0 & Pt198Di==0)
        	
