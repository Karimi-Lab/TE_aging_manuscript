library(stringr)
library(Seurat)


dir <-"/scratch/prj/celgene/CML/J-DU000905/J-DU000905/JGAS000230/JGAD000324/bam_scTE/AddTE"
setwd(dir)
dir2<-"/scratch/prj/celgene/CML/J-DU000905/J-DU000905/JGAS000230/JGAD000324/bam_scTE/MergeTE/"


allMerge <- function(TE) {

  files <- Sys.glob(paste0("*_",TE,".rds"))

  cat(paste0("\033[0;32mReading files\n\033[0m  - ",files[1], "\n"))
  classx <- readRDS(files[1])
  classes <- list()

  for (i in files[-1]){
    writeLines(paste0("  - ",i))
    s <- readRDS(i)
    classes <- c(classes, s)
  }
  all_samples <- merge(classx, classes, merge.data = TRUE)

 # all_samples@meta.data$MajorSep <- as.factor(all_samples@meta.data$MajorSep)
 # all_samples@meta.data$MinorSep <- as.factor(all_samples@meta.data$MinorSep)
  all_samples@meta.data$cell_annotation <- as.factor(all_samples@meta.data$cell_annotation)
  all_samples@meta.data$patient <- as.factor(all_samples@meta.data$patient)
 all_samples@meta.data$group <- as.factor(all_samples@meta.data$group)
  cat("\n\033[0;32mFinding Features\033[0m\n")
  all_samples <- FindVariableFeatures(all_samples)
  all_samples <- ScaleData(all_samples)
  writeLines("\n\033[0;32mRunning PCA\033[0m")
  all_samples <- RunPCA(all_samples, verbose = FALSE)
  writeLines("\n\033[0;32mFinding Neighbours\033[0m")
  all_samples <- FindNeighbors(all_samples, dims = 1:10)
  writeLines("\n\033[0;32mFinding Clusters\033[0m")
  all_samples <- FindClusters(all_samples, resolution = 0.5, verbose = FALSE)
  writeLines("\n\033[0;32mRunning UMAP\033[0m")
  all_samples <- RunUMAP(all_samples, dims = 1:10)

  outfile <- paste0("AllSamples_", TE, ".rds")
  writeLines(c("\n\033[0;32mSaving Seurat object\033[0m",paste0("  - ", outfile)))
  saveRDS(all_samples, paste0(dir2,outfile))
} 

allMerge("class")
allMerge("fam")
