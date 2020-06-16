library(argparse)
library(stringr)
library(flowCore)
library(FlowSOM)
library(Rtsne)
library(ggplot2)
library(sva)
library(harmony)
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--count",
                    type="character",
                    default="output",
                    help="the path to save result")

parser$add_argument("-k",
                    type="integer",
                    default="40")

parser$add_argument("--outdir",
                    type="character",
                    default="")

args <- parser$parse_args()
if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}

data<-readRDS(args$count)
print(paste0("number of cells : ",nrow(data)))
print(tail(rownames(data)))

n_data<-ncol(data)
marker_cols <- 1:n_data 

data_FlowSOM <- flowCore::flowFrame(data)


set.seed(1234)

# run FlowSOM (initial steps prior to meta-clustering)
print("### Run FlowSOM")
pdf(file.path(args$outdir,"FlowSOM.pdf"),width=10,height=12)
out <- FlowSOM::ReadInput(data_FlowSOM, transform = FALSE, scale = FALSE)
out <- FlowSOM::BuildSOM(out, colsToUse = marker_cols)
out <- FlowSOM::BuildMST(out)
dev.off()
# optional visualization

pdf(file.path(args$outdir,"FlowSOM_PlotStars.pdf"),width=10,height=12)
FlowSOM::PlotStars(out)
dev.off()

# extract cluster labels (pre meta-clustering) from output object

labels_pre <- out$map$mapping[, 1]

# specify final number of clusters for meta-clustering (can also be selected 
# automatically, but this often does not perform well)

k <-args$k

# run meta-clustering

# note: In the current version of FlowSOM, the meta-clustering function 
# FlowSOM::metaClustering_consensus() does not pass along the seed argument 
# correctly, so results are not reproducible. We use the internal function 
# ConsensusClusterPlus::ConsensusClusterPlus() to get around this. However, this
# will be fixed in the next update of FlowSOM (version 1.5); then the following 
# (simpler) code can be used instead:
#seed <- 1234
#out <- FlowSOM::metaClustering_consensus(out$map$codes, k = k, seed = seed)
print("### meta Clustering")
seed <- 1234
out <- ConsensusClusterPlus::ConsensusClusterPlus(t(out$map$codes), maxK = k, seed = seed)
out <- out[[k]]$consensusClass


# extract cluster labels from output object

labels <- out[labels_pre]

# summary of cluster sizes and number of clusters

table(labels)
length(table(labels))

# save cluster labels

res <- data.frame(barcode=rownames(data),cluster = labels)

write.table(res, file = file.path(args$outdir,"cluster_labels_FlowSOM.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")



