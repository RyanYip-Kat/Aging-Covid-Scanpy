library(CATALYST)
library(scater)
library(SingleCellExperiment)
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
antigens<-rownames(sce)
for(antigen in antigens){
        jpeg(file.path(args$outdir,paste0(antigen,"_condition_tSNE.jpeg")),width=1024,height=1024)
        p=plotDR(sce,dr="TSNE", color_by =antigen) + facet_wrap("condition")
        print(p)
        dev.off()

        jpeg(file.path(args$outdir,paste0(antigen,"_tSNE.jpeg")),width=1024,height=1024)
        p=plotDR(sce,dr="TSNE", color_by =antigen)
        print(p)
        dev.off()

	jpeg(file.path(args$outdir,paste0(antigen,"_exprs.jpeg")),width=1024,height=1024)
	p=plotExpression(sce,antigen,x = "condition", exprs_values = "exprs",colour = "condition")
        print(p)
	dev.off()
}

