library(argparse)
library(stringr)
library(ggplot2)
library(cicero)

library(data.table)
library(magrittr)

parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--outdir",
                    type="character",
                    default="output",
                    help="the path to save result")

parser$add_argument("--conns",
                    type="character",
                    default="",
                    help="the path of cicero connection from runCicero")

args <- parser$parse_args()

if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}


chr=paste("chr",c(as.character(1:22),"X","Y"))
gene_anno <- rtracklayer::readGFF("Homo_sapiens.GRCh38.99.gtf.gz")
gene_anno$chromosome <- paste0("chr", gene_anno$seqid)
gene_anno$gene <- gene_anno$gene_id
gene_anno$transcript <- gene_anno$transcript_id
gene_anno$symbol <- gene_anno$gene_name

gene_anno=subset(gene_anno,chromosome%in%chr)

conns=read.table(args$conns,header=TRUE)

