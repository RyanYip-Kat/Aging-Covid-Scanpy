library(CATALYST)
library(stringr)
library(argparse)
library(flowCore)
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--fcs",
                    type="character",
                    default="",
                    help="the path of fcs files")


parser$add_argument("--outdir",
                      type="character",
                      default="output",
                      help="save path")

parser$add_argument("--key",
                      type="character",
                      default="output",
                      help="save path")

args<-parser$parse_args()
if(!dir.exists(args$outdir)){
            dir.create(args$outdir,recursive=TRUE)
}

ff<-read.FCS(args$fcs)
name=ff@parameters@data$name
#keys=str_extract(name,"\\d+")
#keys=na.omit(keys)
#keys=as.integer(keys)
#keys<-read.csv("CovID/result/CoVID_Key.csv",stringsAsFactors=F,row.names=1)
keys<-read.csv(args$key,stringsAsFactors=F,row.names=1)
colnames(keys)<-str_extract(colnames(keys),"\\d+")

print("### Assignment of preliminary IDs")
re0 <- assignPrelim(x=ff, y=keys, verbose=FALSE)

print("### Estimation of separation cutoffs")
# estimate separation cutoffs
re <- estCutoffs(x=re0)

print("### Applying deconvolution parameters")
# use population-specific cutoffs
re <- applyCutoffs(x=re)

print("### save fcs")
outFCS(x=re, y=ff,out_path=args$outdir)
saveRDS(re,file.path(args$outdir,"re.rds"))


