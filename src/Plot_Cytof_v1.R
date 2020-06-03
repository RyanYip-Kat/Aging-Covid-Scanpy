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
#pdf(file.path(args$outdir,"Antigen_MedExprs.pdf"),width=12,height=10)
#p <- plotMedExprs(sce, facet = "antigen",
#		  group_by = "condition", shape_by = "patient_id")
#print(p)
#dev.off()

pdf(file.path(args$outdir,"delta_area.pdf"),width=12,height=8)
metadata(sce)$delta_area
dev.off()

pdf(file.path(args$outdir,"plotAbundances_meta16.pdf"),width=12,height=10)
p=plotAbundances(sce, k = "meta16", by = "sample_id", group_by = "condition")
print(p)
dev.off()


#plotExpression
jpeg(file.path(args$outdir,"plotExprs.jpeg"),width=1024,height=1024)
plotExprs(sce)
dev.off()

jpeg(file.path(args$outdir,"clusterHeatmap_abundances.jpeg"),width=1024,height=1024)
plotClusterHeatmap(sce, hm2 = "abundances",
    draw_freqs = TRUE, cluster_anno = FALSE)

dev.off()


