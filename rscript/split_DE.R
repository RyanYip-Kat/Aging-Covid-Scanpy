library(stringr)
library(argparse)

parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--data",
                    type="character",
                    default="",
                    help="")
parser$add_argument("--outdir",
                    type="character",
                    default="output",
                    help="")


args <- parser$parse_args()
if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}

split_status_DE<-function(DATA,outdir){
  if(!dir.exists(outdir)){
    dir.create(outdir,recursive = TRUE)
  }
  stopifnot("cluster"%in%colnames(DATA))
  clusters=unique(DATA$cluster)
  for(c in clusters){
    df1=subset(DATA,cluster==c & pvals<0.05)
    filename1=paste0(c,"_pvals_markers.csv")
    write.table(df1,file.path(outdir,filename1),sep=",",quote = FALSE,row.names = F)
    df2=subset(DATA,cluster==c & pvals_adj<0.05)
    filename2=paste0(c,"_pvals_adj_markers.csv")
    write.table(df2,file.path(outdir,filename2),sep=",",quote = FALSE,row.names = F)

    df3=subset(DATA,cluster==c & logfoldchanges>0.25)
    filename3=paste0(c,"_fc_markers.csv")
    write.table(df3,file.path(outdir,filename3),sep=",",quote = FALSE,row.names = F)

    df4=subset(DATA,cluster==c & logfoldchanges>0.25 & pvals_adj<0.05)
    filename4=paste0(c,"_fc_pvals_adj_markers.csv")
    write.table(df4,file.path(outdir,filename4),sep=",",quote = FALSE,row.names = F)
    
    df5=subset(DATA,cluster==c & logfoldchanges>0.25 & pvals<0.05)
    filename5=paste0(c,"_fc_pvals_markers.csv")
    write.table(df5,file.path(outdir,filename5),sep=",",quote = FALSE,row.names = F)

  }
}

DATA=read.csv(args$data,sep=",")
split_status_DE(DATA,args$outdir)
