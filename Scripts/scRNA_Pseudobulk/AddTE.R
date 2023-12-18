library(stringr)
library(Seurat)

appendTes <- function(asd, sample,TE, meta,l){
  print(paste0(sample,"_",TE))
  mux_class <- readRDS(paste0(dir,sample,"_", TE, "_reads.rds"))
  rownames(mux_class) <- paste0(sapply(rownames(mux_class), function(x) strsplit(x, "-")[[1]][1], USE.NAMES=FALSE),"-",l)
  asd_class <- asd[,colnames(asd) %in% rownames(mux_class)]
  mux_class <- mux_class[rownames(mux_class) %in% colnames(asd_class),]
  meta_class <- meta[rownames(meta) %in% colnames(asd_class),]
  mux_class <- mux_class[match(colnames(asd_class),rownames(mux_class)),]
  asd_class <- rbind(asd_class, t(mux_class))
  so_class <- CreateSeuratObject(asd_class, meta.data = meta_class)
  so_class <- NormalizeData(so_class)
  saveRDS(so_class, paste0(dir,sample, "_", TE, ".rds"))
}

addCellsMeta <- function(CMLSamples){
  cells <- colnames(CMLSamples)
  cells <- cells
  names(cells) <- colnames(CMLSamples)
  CMLSamples <- AddMetaData(CMLSamples, metadata = cells, col.name = "cells")
  return(CMLSamples)
}

CMLSamples <- readRDS("/scratch/prj/celgene/CML/J-DU000905/J-DU000905/JGAS000230/JGAD000324/bam_scTE/AddTE/Japanese_scRNA.rds")
meta <- CMLSamples@meta.data
asd <- GetAssayData(CMLSamples)

dir <-"/scratch/prj/celgene/CML/J-DU000905/J-DU000905/JGAS000230/JGAD000324/bam_scTE/AddTE/"
setwd(dir)
files <- Sys.glob("*_fam_reads.rds")

l=0
for (i in 1:7){
  l=l+1
  print(l)
  print(i)
  sample <- sub("_fam_reads.rds", "", files[i])
  sample2<-str_split(sample,"[.]")[[1]][2]
  print(sample2)
  cell_s<-CMLSamples@meta.data[which(CMLSamples@meta.data$patient==sample2),]$cell_barcode
  asd2 <- asd[,cell_s]
  appendTes(asd2, sample, "class", meta,l)
  appendTes(asd2, sample,"fam", meta,l)
}

l=8
for (i in 8:length(files)){
  l=l+1
  print(l)
  print(i)
  sample <- sub("_fam_reads.rds", "", files[i])
  sample2<-str_split(sample,"[.]")[[1]][2]
  print(sample2)
  cell_s<-CMLSamples@meta.data[which(CMLSamples@meta.data$patient==sample2),]$cell_barcode
  asd2 <- asd[,cell_s]
  appendTes(asd2, sample, "class", meta,l)
  appendTes(asd2, sample,"fam", meta,l)
}



