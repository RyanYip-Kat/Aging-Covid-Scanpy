library(argparse)
library(stringr)
library(flowCore)
library(FlowSOM)
library(Rtsne)
library(ggplot2)
library(sva)

parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--outdir",
                    type="character",
                    default="output",
                    help="the path to save result")


parser$add_argument("--path",
                    type="character",
                    default="")

parser$add_argument("--batch_correct",
                    dest="batch_correct",
                    action="store_true")

parser$add_argument("--pct",
		    type="double",
		    default="1.0")

args <- parser$parse_args()
if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}

cols_file<-file.path(args$path,"cols.csv")
if(file.exists(cols_file)){
	cols<-read.csv(file.path(args$path,"cols.csv"),stringsAsFactors=F)
}else{
	cols<-NULL
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

print("### Loading data")
files<-file.path(args$path,list.files(args$path,pattern=".fcs"))
counts_list<-lapply(files,function(file){
                    prefix<-str_split(basename(file),"\\.")[[1]][1]
		    ff<-flowCore::read.FCS(file,transformation = FALSE, truncate_max_range = FALSE)
		    data <- flowCore::exprs(ff)
		    if(!is.null(cols)){
			    column<-as.character(cols$name)
			    antigen<-as.character(cols$desc)
			    data<-data[,column]
			    colnames(data)<-antigen
		    }
		    row_names<-paste(prefix,1:nrow(data),sep="_")
		    rownames(data)<-row_names
		    if(args$pct!=0){
			    idx<-sample(1:nrow(data),floor(args$pct*nrow(data)),replace=FALSE)
			    data<-data[idx,]
		    }
		    return(data)
		    })
data<-do.call(rbind,counts_list)
saveRDS(data,file.path(args$outdir,paste0(args$pct,"_counts.rds")))
if(args$batch_correct){
	counts<-t(data) # gene x sample
        batch.labels = factor(unlist(lapply(colnames(counts),function(x){return(str_split(x,"_")[[1]][1])})))
        print(table(batch.labels))
	counts_norm = batch.normalise.comBat(counts =counts, batch.groups = batch.labels)
	data<-t(counts_norm)

	tmp=file.path(args$outdir,"batch_corrected")
	if(!dir.exists(tmp)){
          dir.create(tmp,recursive=TRUE)
	}
	saveRDS(data,file.path(tmp,paste0(args$pct,"_counts.rds")))

}

print(paste0("number of cells : ",nrow(data)))
print(tail(rownames(data)))
#file<-args$fcs
#data <- flowCore::exprs(flowCore::read.FCS(file, transformation = FALSE, truncate_max_range = FALSE))

n_data<-ncol(data)
marker_cols <- 1:n_data 

# apply arcsinh transformation
# (with standard scale factor of 5 for CyTOF data; alternatively 150 for flow
# cytometry data; see Bendall et al. 2011, Science, Supplementary Figure S2)

asinh_scale <- 5
data[, marker_cols] <- asinh(data[, marker_cols] / asinh_scale)
# create flowFrame object (required input format for FlowSOM)
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

k <-50

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

res <- data.frame(cluster = labels)

write.table(res, file = file.path(args$outdir,"cluster_labels_FlowSOM.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")




#################
### RUN RTSNE ###
#################

# subsampling (required due to runtime)

n_sub <- nrow(data)

set.seed(1234)
ix <- sample(1:length(labels), n_sub)

# prepare data for Rtsne (matrix format required)

#data_Rtsne <- data[ix, marker_cols]
data_Rtsne <- data[, marker_cols]
data_Rtsne <- as.matrix(data_Rtsne)

dim(data_Rtsne)

# remove any near-duplicate rows (required by Rtsne)

dups <- duplicated(data_Rtsne)
data_Rtsne <- data_Rtsne[!dups, ]

dim(data_Rtsne)
head(rownames(data_Rtsne))

# run Rtsne (Barnes-Hut-SNE algorithm; runtime: 2-3 min)

# note initial PCA is not required, since we do not have too many dimensions
# (i.e. not thousands, which may be the case in other domains)

set.seed(1234)
out_Rtsne <- Rtsne(data_Rtsne, pca = TRUE, verbose = TRUE,num_threads=8)




###################
### CREATE PLOT ###
###################

# load cluster labels (if not still loaded)

file_labels <- file.path(args$outdir,"cluster_labels_FlowSOM.txt")
data_labels <- read.table(file_labels, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
labels <- data_labels[, "cluster"]

# select points used by Rtsne

#labels_plot <- labels[ix][!dups]
labels_plot <- labels[!dups]
length(labels_plot)  ## should be same as number of rows in data_Rtsne

# prepare Rtsne output data for plot

data_plot <- as.data.frame(out_Rtsne$Y)
colnames(data_plot) <- c("tSNE_1", "tSNE_2")

head(data_plot)
dim(data_plot)  ## should match length of labels_plot (otherwise labels will not match up correctly)

data_plot[, "cluster"] <- as.factor(labels_plot)
rownames(data_plot)<-rownames(data_Rtsne)
status=unlist(lapply(rownames(data_plot),function(x){return(str_split(x,"_")[[1]][1])}))
data_plot[,"status"]<-as.factor(status)
head(data_plot)

saveRDS(data_plot,file.path(args$outdir,"tSNE.rds"))
# plot 2-dimensional t-SNE projection
pdf(file.path(args$outdir,"FlowSOM_Rtsne_Clusters.pdf"),width=12,height=8)
ggplot(data_plot, aes(x = tSNE_1, y = tSNE_2, color = cluster)) + 
  geom_point(size = 0.2) + 
  coord_fixed(ratio = 1) + 
  ggtitle("t-SNE projection with FlowSOM clustering") + 
  theme_bw()

dev.off()

pdf(file.path(args$outdir,"FlowSOM_Rtsne_Status.pdf"),width=12,height=8)
ggplot(data_plot, aes(x = tSNE_1, y = tSNE_2, color = status)) +
    geom_point(size = 0.2) +
    coord_fixed(ratio = 1) +
    ggtitle("t-SNE projection with FlowSOM clustering") +
    theme_bw()

dev.off()

