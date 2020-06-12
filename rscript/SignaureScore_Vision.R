library(argparse)
library(stringr)
library(Seurat)
library(VISION)

#############################
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--seurat",
                    type="character",
                    default="")

parser$add_argument("--column",
                    type="character",
                    default=NULL)

parser$add_argument("--subset",
		    nargs="+",
                    type="character",
                    default=NULL)

parser$add_argument("--projection",
                    type="character",
                    default=NULL)

parser$add_argument("--outdir",
                    type="character",
                    default="vision_results")

parser$add_argument("--sig",
                    type="character",
                    default="data/h.all.v7.1.symbols.gmt")

args <- parser$parse_args()

if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}

options(mc.cores = 4)

seurat_obj<-readRDS(args$seurat)
DefaultAssay(seurat_obj)="RNA"
metadata=seurat_obj@meta.data
meta=metadata[,"ident",drop=F]
rownames(meta)=rownames(metadata)
if(!is.null(args$column)){
	if(!args$column%in%colnames(metadata)){
		stop("Invaild columns in metadata!")
	}else{
		s<-paste(args$subset,collapse=",")
		print(paste0("Subset :",args$column," with :",s))
		Idents(seurat_obj)=metadata[[args$column]]
		seurat_obj<-subset(seurat_obj,idents=args$subset)
	}
}

print(paste0("Using signatures from :",args$sig))
signatures <-args$sig

print(" Analysis ")
n_cells=ncol(seurat_obj)
seurat_reduction=names(seurat_obj@reductions)
if(!is.null(args$projection)){
	if(!args$projection%in%seurat_reduction){
		stop("Invaild projection input!")
	}
	vision.obj <- Vision(seurat_obj,
			     signatures = signatures,
			     meta=meta,
			     dimRed=args$projection,
                             projection_methods =NULL)
}else{
	vision.obj <- Vision(seurat_obj,
			     meta=meta,
                             signatures = signatures)
}
vision.obj <- analyze(vision.obj)
saveRDS(vision.obj, file.path(args$outdir,'vision_results.rds'))

#print("View result")
#viewResults(vision.obj,host="10.100.110.103",port="1234",browser=FALSE)



