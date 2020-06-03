library(argparse)
library(harmony)
library(CATALYST)
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

parser$add_argument("--number",
                    type="integer",
                    default="20000")
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

#if(!is.null(cols)){
#        column<-as.character(cols$name)
#        column=ifelse(str_detect(column,"\\d+"),str_extract(column,"\\d+"),column)
#        antigen<-as.character(cols$marker)
#        data<-data[,column]
#        colnames(data)<-antigen
#}


print("### Loading data")
files<-list.files(args$path,pattern=".fcs",recursive=TRUE,full.names=TRUE)
print(files)
counts_list<-lapply(files,function(file){
	            p1<-basename(dirname(file))
                    p2<-str_split(basename(file),"\\.")[[1]][1]
		    prefix<-paste(p1,p2,sep="-")
		    ff<-flowCore::read.FCS(file,transformation = FALSE, truncate_max_range = FALSE)
		    data <- flowCore::exprs(ff)
		    Names=colnames(data)
		    col_names=ifelse(str_detect(Names,"\\d+"),str_extract(Names,"\\d+"),Names)
		    colnames(data)=col_names
		    row_names<-paste(prefix,1:nrow(data),sep="_")
		    rownames(data)<-row_names
		    if(!is.null(cols)){
		    	    column<-as.character(cols$name)
                            column=ifelse(str_detect(column,"\\d+"),str_extract(column,"\\d+"),column)
                            antigen<-as.character(cols$marker)
                            data<-data[,column]
                            colnames(data)<-antigen
		    }

		    if(args$number!=0){
			    idx<-sample(1:nrow(data),as.integer(args$number),replace=FALSE)
			    data<-data[idx,]
		    }
		    print(dim(data))
		    return(data)
		    })

print("Merge data")
data<-do.call(rbind,counts_list)
print(paste0("number of cells : ",nrow(data)))
print(tail(rownames(data)))
print(dim(data))

Names<-rownames(data)
#file<-args$fcs
#data <- flowCore::exprs(flowCore::read.FCS(file, transformation = FALSE, truncate_max_range = FALSE))
#data_FlowSOM <- flowCore::flowFrame(data)
#print("### CATALYST : normCytof")
#pdf(file.path(args$outdir,"normCytof.pdf"),width=12,height=10)
#CATALYST::normCytof(x=data_FlowSOM,y="dvs",remove_beads=FALSE,out_path=args$outdir,plot=TRUE)
#dev.off()


#file<-list.files(getwd(),"normalised.fcs")
#if(length(file)!=0){
#	system(paste0("mv ",file," ",args$outdir))
#}

#normalised.fcs<-file.path(args$outdir,list.files(args$outdir,"normalised.fcs"))
#data_FlowSOM=flowCore::read.FCS(normalised.fcs)
#data<-flowCore::exprs(data_FlowSOM)

#if(!is.null(cols)){
#	column<-as.character(cols$name)
#	column=ifelse(str_detect(column,"\\d+"),str_extract(column,"\\d+"),column)
#        antigen<-as.character(cols$marker)
#        data<-data[,column]
#        colnames(data)<-antigen
#}
#print(dim(data))
#rownames(data)<-Names
saveRDS(data,file.path(args$outdir,"counts.rds"))
status=unlist(lapply(Names,function(x){return(str_split(x,"-")[[1]][1])}))
sample=unlist(lapply(Names,function(name){return(str_split(name,"_")[[1]][1])}))
print(table(sample))
print(table(status))
if(args$batch_correct){
        metadata<-data.frame("cell_id"=Names,"sample"=sample,"status"=status)
        harmony_embeddings<-HarmonyMatrix(as.matrix(data), metadata, c("sample")) #v3
        rownames(harmony_embeddings)<-Names
        colnames(harmony_embeddings)=paste("Harmony",1:ncol(harmony_embeddings),sep="_")
        saveRDS(harmony_embeddings,file.path(args$outdir,"harmony.rds"))
}


n_data<-ncol(data)
marker_cols <- 1:n_data 

# apply arcsinh transformation
# (with standard scale factor of 5 for CyTOF data; alternatively 150 for flow
# cytometry data; see Bendall et al. 2011, Science, Supplementary Figure S2)

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

k <-30

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
print("### Run TSNE")
p1=proc.time()
out_Rtsne <- Rtsne(data_Rtsne, pca = FALSE, verbose = TRUE,num_threads=4,max_iter=2000)
p2=proc.time()
time_delta<-p2[3]-p1[3]
sprintf("run TSNE use %f seconds",time_delta)



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

head(data_plot)
status=unlist(lapply(rownames(data_plot),function(x){return(str_split(x,"-")[[1]][1])}))
data_plot[,"status"]<-as.factor(status)
saveRDS(data_plot,file.path(args$outdir,"tSNE.rds"))
# plot 2-dimensional t-SNE projection
pdf(file.path(args$outdir,"FlowSOM_Rtsne_Cluster.pdf"),width=12,height=8)
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
