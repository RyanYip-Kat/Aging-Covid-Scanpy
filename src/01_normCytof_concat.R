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

fcs_path<-file.path(args$path,list.files(args$path,pattern=".fcs",recursive=TRUE))
process<-function(fcs){
         ff<-flowCore::read.FCS(fcs)
         pdf(file.path(args$outdir,paste0(basename(fcs),".pdf")),width=10,height=8)
         normCytof(x=ff, y="dvs", 
                   out_path=args$outdir,
		   k=120, plot=TRUE)
	 dev.off()
}


ff <- concatFCS(fcs_path,file_num =FALSE)
#DATA=ff@exprs
#filename=DATA[,"FileNum"]
#cells=paste("S",1:nrow(DATA),sep="-")
#cells=paste(cells,filename,sep="_")
#rownames(DATA)=cells
#DATA<-DATA[,-ncol(DATA)]
#ff@exprs=DATA
#data<-ff@parameters@data
#ff@parameters@data=data[-nrow(data),]
#ff <- concatFCS(fcs_path)
pdf(file.path(args$outdir,"concat.pdf"),width=10,height=8)
normCytof(x=ff, y="dvs",
          out_path=args$outdir,
          plot=TRUE)
dev.off()

