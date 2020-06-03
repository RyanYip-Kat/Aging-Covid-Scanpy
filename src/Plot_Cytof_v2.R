library(flowCore)
library(CATALYST)
library(scater)
library(SingleCellExperiment)
library(Seurat)
library(ggplot2)
library(argparse)

parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--outdir",
                    type="character",
                    default="output",
                    help="the path to save result")


parser$add_argument("--sce",
		    type="character",
		    default="")
args <- parser$parse_args()
if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}

sce<-readRDS(args$sce)
antigens=rownames(sce)
for(antigen in antigens){
	p=plotDR(sce,dr="TSNE", color_by =antigen)+ facet_wrap("patient_id")
	pdf(file.path(args$outdir,paste0(antigen,"_tSNE_","patient_id.pdf")),width=12,height=10)
	print(p)
	dev.off()
}
