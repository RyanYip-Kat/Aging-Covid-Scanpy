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
### Compute TPM expression values from raw UMI counts
tpm <- function(counts, mult=10000,eps=1e-6)
{
  info("Running TPM normalisation")
  total.counts = colSums(counts)
  scaled.counts = t(t(counts) / total.counts)
  scaled.counts * mult
}
### Run ComBat batch correction from the SVA package
batch.normalise.comBat <- function(counts, batch.groups, max.val=6)
{
  batch.groups = factor(batch.groups) ## drop zero levels
  batch.id = 1:length(unique(batch.groups))
  names(batch.id) = unique(batch.groups)
  batch.ids = batch.id[batch.groups]
  correct.data = ComBat(counts,batch.ids, prior.plots=FALSE, par.prior=T)
  correct.data[correct.data > max.val] = max.val
  as.data.frame(correct.data)
}

counts<-t(counts) # gene x sample
batch.labels = factor(unlist(lapply(colnames(counts),function(x){return(str_split(x,"_")[[1]][1])})))
print(table(batch.labels))

counts_norm = batch.normalise.comBat(counts =counts, batch.groups = batch.labels)

saveRDS(counts_norm,file.path(args$outdir,"counts_norm.rds"))
