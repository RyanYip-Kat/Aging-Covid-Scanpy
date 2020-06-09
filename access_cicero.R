library(cicero)
cicero_conn.file = paste0(down_dir, '/cicero_interactions.txt')
Cicero_Plot_Region = 'chr5:140610000-140640000'
if(file.exists(cicero_conn.file)){
  conns = fread(cicero_conn.file)
  conns = data.frame(conns)
  temp <- tempfile()
  if(grepl(GENOME_NAME, pattern = 'mm10', ignore.case = T)) {
    download.file('ftp://ftp.ensembl.org/pub/release-95/gtf/mus_musculus/Mus_musculus.GRCm38.95.gtf.gz', temp)
  }
  
  if(grepl(GENOME_NAME, pattern = 'mm9', ignore.case = T)) {
    download.file('ftp://ftp.ensembl.org/pub/release-67/gtf/mus_musculus/Mus_musculus.NCBIM37.67.gtf.gz', temp)
  }
  
  if(grepl(GENOME_NAME, pattern = 'hg38', ignore.case = T)) {
    download.file('ftp://ftp.ensembl.org/pub/release-95/gtf/homo_sapiens/Homo_sapiens.GRCh38.95.gtf.gz', temp)
  }
  
  if(grepl(GENOME_NAME, pattern = 'hg19', ignore.case = T)) {
    download.file('ftp://ftp.ensembl.org/pub/release-67/gtf/homo_sapiens/Homo_sapiens.GRCh37.67.gtf.gz', temp)
  }
  
  gene_anno <- rtracklayer::readGFF(temp)
  unlink(temp)
  
  # rename some columns to match requirements
  gene_anno$chromosome <- paste0("chr", gene_anno$seqid)
  gene_anno$gene <- gene_anno$gene_id
  gene_anno$transcript <- gene_anno$transcript_id
  gene_anno$symbol <- gene_anno$gene_name
  gene_anno = subset(gene_anno, select = c(chromosome, start, end, strand, 
                                           transcript, gene, symbol))
  gene_anno = gene_anno[complete.cases(gene_anno), ]
  
  chr0 = unlist(strsplit(Cicero_Plot_Region, ':'))[1] ## chr5:140610000-140640000
  region0 = unlist(strsplit(Cicero_Plot_Region, ':'))[2]
  start0 = as.integer(unlist(strsplit(region0, '-'))[1])
  end0 = as.integer(unlist(strsplit(region0, '-'))[2])
  
  
  cicero::plot_connections(conns, chr0, start0, end0,
                   gene_model = gene_anno, 
                   coaccess_cutoff = .3, 
                   connection_width = 1, 
                   collapseTranscripts = "longest",
                   viewpoint_alpha = 0)
   
}
