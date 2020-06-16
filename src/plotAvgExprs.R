library(ggplot2)
library(Seurat)
library(reshape2)
library(argparse)

parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--outdir",
                    type="character",
                    default="output",
                    help="the path to save result")


parser$add_argument("--seurat",
                    type="character",
                    default="")
args <- parser$parse_args()
if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}

seurat<-readRDS(args$seurat)
markers=rownames(seurat)
Idents(seurat)=seurat$condition
conditions=unique(seurat$condition)

thm <-  theme(axis.text.x = element_text(
            angle = 45, hjust = 1, vjust = 1))

thm<-thm+ylab("expression") + theme_bw() + thm + theme(
            legend.key.height  =  unit(0.8, "lines"),
            axis.text = element_text(color = "black"),
            strip.text = element_text(face = "bold"),
            strip.background = element_rect(fill = NA, color = NA),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_line(color = "grey", size = 0.2))


plot_list<-lapply(1:length(markers),function(k){
	       marker<-markers[k]
	       avg_list<-list()
	       for(i in 1:length(conditions)){
		       x<-subset(seurat,idents=conditions[i])
		       Idents(x)<-x$sample_id
		       avg=AverageExpression(x,features=marker,verbose=FALSE)[["exprs"]]
		       avg$condition=conditions[i]
		       avg=melt(avg)
		       avg_list[[i]]=avg
	       }
	       mat<-do.call(rbind,avg_list)
	       p=ggplot(mat,aes(condition,value,col=condition))+
		       geom_point(alpha = 0.8,aes(fill =condition, shape =variable),
				  position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0))+
                       geom_boxplot(alpha = 0.4, width = 0.8,fill = NA, outlier.color = NA, show.legend = FALSE)+                       ggtitle(marker)+thm
	       
	       pdf(file.path(args$outdir,paste0(marker,".pdf")),width=12,height=10)
	       print(p)
	       dev.off()
	     })

#pdf(file.path(args$outdir,"ExprAntigen.pdf"),width=12,height=10)
#print(CombinePlots(plot_list,ncol=4))
#dev.off()

