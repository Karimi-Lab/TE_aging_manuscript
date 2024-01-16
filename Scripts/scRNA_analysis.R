##### Setup #####
`%notin%` <- Negate(`%in%`)

# BiocManager::install(c("apeglm","clusterProfiler","DESeq2","vsn"))
library(AnnotationDbi)
library(ape)
library(apeglm)
library(circlize)
library(clusterProfiler)
library(ComplexHeatmap)
library(cowplot)
library(data.table)
library(DESeq2)
library(dplyr)
library(edgeR)
library(export)
library(finalfit)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(ggstatsplot)
library(ggthemes)
library(gmodels)
library(grid)
library(GSVA)
library(limma)
library(MASS)
library(Matrix)
library(msigdbr)
library(NMF)
library(nnet)
library(org)
library(org.Hs.eg.db)
library(patchwork)
library(phangorn)
library(pheatmap)
library(psych)
library(RColorBrewer)
library(readr)
library(readxl)
library(reshape2)
library(scales)
library(scatterplot3d)
library(Seurat)
library(singscore)
library(stringi)
library(stringr)
library(survival)
library(survminer)
library(tidyverse)
library(vcd)
library(vsn)


##scRNA-seq #####
setwd("../Data/Single_Cell/scRNASeq-PBMC/")
pseudo_class<-readRDS("./pseudo_class.rds")
pseudo_fam<-readRDS("./pseudo_fam.Rds")

##class
unique(pseudo_class$Sample_ID)
young_sample<-c("A06", "A08", "A10", "A11" ,"A12", "A18" ,"A20" ,"A21" ,"A23" ,"A24" ,"A26")
old_sample<-c("D03", "D15", "D22" ,"D24", "E04" ,"E05", "E08" ,"E10", "E16" ,"E17")
pseudo_class$Status<-NA
pseudo_class$Status<-"young"
pseudo_class@meta.data[which(pseudo_class@meta.data$Sample_ID%in%old_sample),]$Status<-"old"
saveRDS(pseudo_class,"./pseudo_class.rds")
##fam
unique(pseudo_fam$Sample_ID)
pseudo_fam$Status<-NA
pseudo_fam$Status<-"young"
pseudo_fam@meta.data[which(pseudo_fam@meta.data$Sample_ID%in%old_sample),]$Status<-"old"
saveRDS(pseudo_fam,"./pseudo_fam.rds")

##fam
unique(pseudo_fam$Sample_ID)
A<-pseudo_fam@meta.data
#write.table(A,file = "./pseudo_fam_metadata.txt",sep="\t",quote = F)
pseudo_fam_metadata <- read_delim("pseudo_fam_metadata.txt", 
                                    delim = "\t", escape_double = FALSE, 
                                    trim_ws = TRUE)
pseudo_fam_metadata<-as.data.frame(pseudo_fam_metadata)
rownames(pseudo_fam_metadata)<-pseudo_fam_metadata$...1
pseudo_fam_metadata<-pseudo_fam_metadata[,-1]
pseudo_fam@meta.data<-pseudo_fam_metadata
pseudo_fam$Status<-NA
pseudo_fam$Status<-"young"
pseudo_fam@meta.data[which(pseudo_fam@meta.data$Sample_ID%in%old_sample),]$Status<-"old"
# saveRDS(pseudo_fam,"./pseudo_fam.rds")

############### scRNA-seq MUX plot
setwd("../Data/Single_Cell/scRNASeq-PBMC/")
pseudo_class<-readRDS("./pseudo_class.rds")
pseudo_fam<-readRDS("./pseudo_fam.Rds")
TE_masker <- read.csv("../Data/Single_Cell/hg38_repeat_masker.gz", header=TRUE, sep="\t")
TE_masker <- unique(TE_masker[,c("repName","repClass","repFamily")])
Class_for_comparison<-unique(TE_masker$repClass)
fam_for_comparison<-unique(TE_masker$repFamily)
pseudo_class_matrix<-as.matrix(pseudo_class@assays$RNA@counts)
pseudo_class_metadata<-pseudo_class@meta.data
pseudo_class_metadata<-pseudo_class_metadata[colnames(pseudo_class_matrix),]
pseudo_class_matrix_normalized<-pseudo_class_matrix[rownames(pseudo_class_matrix)[which(rownames(pseudo_class_matrix)%in%Class_for_comparison)],]
pseudo_class_matrix_normalized_done<-matrix(data = NA,nrow = nrow(pseudo_class_matrix_normalized),ncol = ncol(pseudo_class_matrix_normalized))
rownames(pseudo_class_matrix_normalized_done)<-rownames(pseudo_class_matrix_normalized)
colnames(pseudo_class_matrix_normalized_done)<-colnames(pseudo_class_matrix_normalized)
#View(pseudo_class_matrix_normalized)
#delete<-data.frame(sum=colSums(pseudo_class_matrix_normalized))

for (i in 1:ncol(pseudo_class_matrix_normalized)) {
  pseudo_class_matrix_normalized_done[,i]<-
    1000000*pseudo_class_matrix_normalized[,i]/(pseudo_class_metadata[i,]$nCount_RNA)
                                                                               #-as.numeric(delete[i,]))
}

pseudo_class_matrix_normalized_done<-t(pseudo_class_matrix_normalized_done)
pseudo_class_metadata<-pseudo_class_metadata[rownames(pseudo_class_matrix_normalized_done),]
pseudo_class_metadata<-cbind(pseudo_class_metadata,pseudo_class_matrix_normalized_done)
pseudo_class_metadata$Status<-factor(pseudo_class_metadata$Status,levels = c("young","old"))
Cell_type<-unique(pseudo_class_metadata$Cell_type)
class<-unique(rownames(pseudo_class_matrix_normalized))

for (i in 1:length(Cell_type)) {
  for (l in 1:length(class)) {
tmp<-pseudo_class_metadata[which(pseudo_class_metadata$Cell_type==Cell_type[i]),]
ggplot(tmp,aes(x=Status,y=log10(tmp[,class[l]]),color=Status))+
  stat_boxplot(geom = "errorbar",width=0.2,aes(color=Status))+
  geom_boxplot(notch = F,width=0.6,outlier.shape = NA)+
  scale_color_manual(values=c("#999999","#1C769F"))+
  labs(x = NULL,y = paste0(Cell_type[i],": log10(RPM) class: ",class[l]))+
  stat_summary(fun="mean",color="black")+
  geom_jitter(shape=16,position = position_jitter(0.1),color="black",size=3)+
  stat_compare_means()+
  theme_classic2()+
  theme(axis.title.x = element_blank(),
  axis.ticks.x=element_blank(),
  axis.title.y = element_text(size = 16, face = "bold"),
  axis.text.y = element_text(size = 16),
   legend.text = element_blank(),
   legend.title = element_blank())
ggsave(paste0("./plot/class/",class[l],"/",as.character(Cell_type[i]),".png"),width = 6,height=8)
}
}

### fam
pseudo_fam_matrix<-as.matrix(pseudo_fam@assays$RNA@counts)
pseudo_fam_metadata<-pseudo_fam@meta.data
pseudo_fam_metadata<-pseudo_fam_metadata[colnames(pseudo_fam_matrix),]
pseudo_fam_matrix_normalized<-pseudo_fam_matrix[rownames(pseudo_fam_matrix)[which(rownames(pseudo_fam_matrix)%in%fam_for_comparison)],]
pseudo_fam_matrix_normalized_done<-matrix(data = NA,nrow = nrow(pseudo_fam_matrix_normalized),ncol = ncol(pseudo_fam_matrix_normalized))
rownames(pseudo_fam_matrix_normalized_done)<-rownames(pseudo_fam_matrix_normalized)
colnames(pseudo_fam_matrix_normalized_done)<-colnames(pseudo_fam_matrix_normalized)
#delete<-data.frame(sum=colSums(pseudo_fam_matrix_normalized))

for (i in 1:ncol(pseudo_fam_matrix_normalized)) {
  pseudo_fam_matrix_normalized_done[,i]<-1000000*pseudo_fam_matrix_normalized[,i]/(pseudo_fam_metadata[i,]$nCount_RNA)
                                                                           #-as.numeric(delete[i,]))
}
pseudo_fam_matrix_normalized_done<-t(pseudo_fam_matrix_normalized_done)
pseudo_fam_metadata<-pseudo_fam_metadata[rownames(pseudo_fam_matrix_normalized_done),]
pseudo_fam_metadata<-cbind(pseudo_fam_metadata,pseudo_fam_matrix_normalized_done)
pseudo_fam_metadata$Status<-factor(pseudo_fam_metadata$Status,levels = c("young","old"))
Cell_type<-unique(pseudo_fam_metadata$Cell_type)
fam<-unique(rownames(pseudo_fam_matrix_normalized))

for (i in 1:length(Cell_type)) {
  for (l in 1:length(fam)) {
    tmp<-pseudo_fam_metadata[which(pseudo_fam_metadata$Cell_type==Cell_type[i]),]
    ggplot(tmp,aes(x=Status,y=log10(tmp[,fam[l]]),color=Status))+
      stat_boxplot(geom = "errorbar",width=0.2,aes(color=Status))+
      geom_boxplot(notch = F,width=0.6,outlier.shape = NA)+
      scale_color_manual(values=c("#999999","#1C769F"))+
      labs(x = NULL,y = paste0(Cell_type[i],": log10(RPM) fam: ",fam[l]))+
      stat_summary(fun="mean",color="black")+
      geom_jitter(shape=16, position = position_jitter(0.1),color="black",size=3)+
      stat_compare_means()+theme_classic2()+
      theme(axis.title.x = element_blank(),
            axis.ticks.x=element_blank(),
            axis.title.y = element_text(size = 16, face = "bold"),
            axis.text.y = element_text(size = 16),
            legend.text = element_blank(),
            legend.title = element_blank())
    ggsave(filename =paste0("./plot/fam/",fam[l],"/",as.character(Cell_type[i]),".png"),width = 6,height=8)
  }
}

##### summary plot
my_Plot_class <- function(cell_type){
  tmp<-pseudo_class_metadata[which(pseudo_class_metadata$Cell_type==cell_type),]
  ggplot(tmp,aes(x=Status,y=log10(tmp[,classes]),color=Status))+
    stat_boxplot(geom = "errorbar",width=0.2,aes(color=Status))+
    geom_boxplot(notch = F,width=0.6,outlier.shape = NA)+
    scale_color_manual(values=c("#999999","#1C769F"))+
    labs(x = NULL,
         title = paste0(cell_type),
         #y =paste0("log10(RPM)")
         y =NULL)+
    stat_summary(fun="mean",color="black")+
    geom_jitter(shape=16, position = position_jitter(0.1),color="black",size=3)+
    stat_compare_means()+theme_classic2()+
    theme(axis.title.x = element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.x=element_blank(),
          axis.title.y = element_text(size = 12, face = "bold"),
          axis.text.y = element_text(size = 12),
          legend.text = element_blank(),
          legend.title = element_blank(),
          legend.position = "none",
          plot.title = element_text(size = 14, color = "black",face = "bold"))
}

my_Plot_fam <- function(cell_type){
  tmp<-pseudo_fam_metadata[which(pseudo_fam_metadata$Cell_type==cell_type),]
  ggplot(tmp,aes(x=Status,y=log10(tmp[,family]),color=Status))+
    stat_boxplot(geom = "errorbar",width=0.2,aes(color=Status))+
    geom_boxplot(notch = F,width=0.6,outlier.shape = NA)+
    scale_color_manual(values=c("#999999","#1C769F"))+
    labs(x = NULL,
         title = paste0(cell_type),
         #y =paste0("log10(RPM)")
         y =NULL)+
    stat_summary(fun="mean",color="black")+
    geom_jitter(shape=16, position = position_jitter(0.1),color="black",size=3)+
    stat_compare_means()+theme_classic2()+
    theme(axis.title.x = element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.x=element_blank(),
          axis.title.y = element_text(size = 12, face = "bold"),
          axis.text.y = element_text(size = 12),
          legend.text = element_blank(),
          legend.title = element_blank(),
          legend.position = "none",
          plot.title = element_text(size = 14, color = "black",face = "bold"))
}

for (l in 1:length(class)) {
  classes<-class[l]
  pdf(file =paste0("./plot/class/Summarised_",classes,".pdf"),width = 15,height=20)
  my_Plot_List_class <- lapply(Cell_type,my_Plot_class)
  pt<-plot_grid(plotlist = my_Plot_List_class, align = "hv", 
            nrow = 5,ncol = 5)
  print(pt)
  dev.off()
}

for (l in 1:length(fam)) {
  family<-fam[l]
  pdf(file =paste0("./plot/fam/Summarised_",family,".pdf"),width = 15,height=20)
  my_Plot_List_fam <- lapply(Cell_type,my_Plot_fam)
  pt<-plot_grid(plotlist = my_Plot_List_fam, align = "hv", 
            nrow = 5,ncol = 5)
  print(pt)
  dev.off()
}

########################### Inflammation ########################
setwd("../Data/Single_Cell/Inflammatory_analysis/")
pseudo_scRNA<-readRDS("../Data/Single_Cell/scRNASeq-PBMC/pseudo_class.rds")
#pseudo_scATAC<-readRDS("/Users/zouxinchen/Desktop/KCL/Minor_project/CML/scATAC-seq_PBMC/pseudo_class.rds")
TE_masker <- read.csv("../Data/Single_Cell/hg38_repeat_masker.gz", header=TRUE, sep="\t")
uniq_TE_masker <- unique(TE_masker[,c("repName","repClass","repFamily")])
Class_for_comparison <- unique(uniq_TE_masker$repClass)


##### inflammation genelists #####
Inflamm_gene_list <- read_excel("../Data/Single_Cell/Inflammatory_analysis/Supplementary_tables (1).xlsx",sheet = "S1")
Inflamm_gene_list <- as.data.frame(Inflamm_gene_list)
Lists<-colnames(Inflamm_gene_list)
subset.matrix <- as.matrix(pseudo_scRNA@assays$RNA@counts)[which(rownames(as.matrix(pseudo_scRNA@assays$RNA@counts))%notin%c("LINE","SINE","LTR")), ]

pseudo_class <- readRDS("../Data/Single_Cell/scRNASeq-PBMC/pseudo_class.rds")
pseudo_fam <- readRDS("../Data/Single_Cell/Xinchen/scRNASeq-PBMC/pseudo_fam.Rds")
for (z in 1:length(Lists)) {
  # z=1
  geneset_used<-Lists[z]
  print(Lists[z])
  Inflamm_genes<-unique(as.character(Inflamm_gene_list[,geneset_used]))[which(!is.na(unique(as.character(Inflamm_gene_list[,geneset_used]))))]
  # Pull the raw expression matrix from the original Seurat object containing only the genes of interest
  object2 <- CreateSeuratObject(subset.matrix) # Create a new Seurat object with just the genes of interest
  orig.ident <- pseudo_scRNA@meta.data # Pull the identities from the original Seurat object as a data.frame
  object2 <- AddMetaData(object = object2, metadata = orig.ident) # Add the idents to the meta.data slot
   # Assign identities for the new Seurat object
  Delete<-colSums(as.matrix(pseudo_scRNA@assays$RNA@counts)[c("SINE","LINE","LTR"),])
  object2$nCount_RNA<-object2$nCount_RNA-Delete
  Matrix_patient<-matrix(NA,nrow = nrow(subset.matrix),ncol = length(unique(object2@meta.data$Sample_ID)))
  rownames(Matrix_patient)<-rownames(subset.matrix)
  colnames(Matrix_patient)<-unique(object2@meta.data$Sample_ID)
  Patient_meta<-data.frame(sample=unique(object2@meta.data$Sample_ID),ncount=NA)
  
  for (i in unique(object2@meta.data$Sample_ID)) {
    ident<-rownames(object2@meta.data[which(object2@meta.data$Sample_ID==i),])
    Matrix_patient[,i]<-rowSums(subset.matrix[,ident])
    Patient_meta[which(Patient_meta$sample==i),]$ncount<-sum(subset.matrix[,ident])
    }
  
  Matrix_patient_normalized_done<-
    matrix(data = NA,nrow = nrow(Matrix_patient),ncol = ncol(Matrix_patient))
  rownames(Matrix_patient_normalized_done)<-rownames(Matrix_patient)
  colnames(Matrix_patient_normalized_done)<-colnames(Matrix_patient)
  
  for (i in 1:ncol(Matrix_patient)) {
    Matrix_patient_normalized_done[,i]<-
      1000000*Matrix_patient[,i]/(Patient_meta[i,]$ncount)
  }
  
  rankData<-rankGenes(Matrix_patient_normalized_done)
  tgfb_gs_up<-as.character(Inflamm_genes)
  scoredf <- simpleScore(rankData, upSet = tgfb_gs_up,knownDirection = T,centerScore = F,dispersionFun = mad)
  High_inflammatory<-rownames(scoredf[which(scoredf$TotalScore>=median(scoredf$TotalScore)),])
  Low_inflammatory<-rownames(scoredf[which(scoredf$TotalScore<median(scoredf$TotalScore)),])
  
  paste0("High: ",High_inflammatory)
  paste0("Low: ",Low_inflammatory)

  Class_for_comparison <- unique(uniq_TE_masker$repClass)
  fam_for_comparison <- unique(uniq_TE_masker$repFamily)
  pseudo_class_matrix<-as.matrix(pseudo_class@assays$RNA@counts)
  pseudo_class_metadata<-pseudo_class@meta.data
  pseudo_class_metadata<-pseudo_class_metadata[colnames(pseudo_class_matrix),]
  pseudo_class_matrix_normalized<-pseudo_class_matrix[rownames(pseudo_class_matrix)[which(rownames(pseudo_class_matrix)%in%Class_for_comparison)],]
  pseudo_class_matrix_normalized_done<-matrix(data = NA,nrow = nrow(pseudo_class_matrix_normalized),ncol = ncol(pseudo_class_matrix_normalized))
  rownames(pseudo_class_matrix_normalized_done)<-rownames(pseudo_class_matrix_normalized)
  colnames(pseudo_class_matrix_normalized_done)<-colnames(pseudo_class_matrix_normalized)
  #View(pseudo_class_matrix_normalized)
  #delete<-data.frame(sum=colSums(pseudo_class_matrix_normalized))
  
  for (i in 1:ncol(pseudo_class_matrix_normalized)) {
    pseudo_class_matrix_normalized_done[,i]<-
      1000000*pseudo_class_matrix_normalized[,i]/(pseudo_class_metadata[i,]$nCount_RNA)
    #-as.numeric(delete[i,]))
  }
  pseudo_class_matrix_normalized_done<-t(pseudo_class_matrix_normalized_done)
  pseudo_class_metadata<-pseudo_class_metadata[rownames(pseudo_class_matrix_normalized_done),]
  pseudo_class_metadata<-cbind(pseudo_class_metadata,pseudo_class_matrix_normalized_done)
  pseudo_class_metadata$Status<-factor(pseudo_class_metadata$Status,levels = c("young","old"))
  pseudo_class_metadata$Inflammatory_Status<-"low_inflammation"
  pseudo_class_metadata[which(pseudo_class_metadata$Sample_ID%in%High_inflammatory),]$Inflammatory_Status<-"high_inflammation"
  pseudo_class_metadata$Inflammatory_Status<-factor(pseudo_class_metadata$Inflammatory_Status,levels = c("low_inflammation","high_inflammation"))
  Cell_type<-unique(pseudo_class_metadata$Cell_type)
  class<-unique(rownames(pseudo_class_matrix_normalized))
  # No text plots
  for (i in 1:length(Cell_type)) {
    for (l in 1:length(class)) {
      tmp<-pseudo_class_metadata[which(pseudo_class_metadata$Cell_type==Cell_type[i]),]
      pdf(file = paste0("./plot_notext/class/",class[l],"/",geneset_used,"/",as.character(Cell_type[i]),".pdf"),
          width = 88*0.0394/3,
          height = 88*0.0394/3)
      pt <- ggplot(tmp,aes(x=Inflammatory_Status,y=log10(tmp[,class[l]]),color=Inflammatory_Status))+
        stat_boxplot(geom = "errorbar",width=0.2,aes(color=Inflammatory_Status))+
        geom_boxplot(notch = F,width=0.6,outlier.shape = NA)+
        scale_color_manual(values=c("#999999","#1C769F"))+
        labs(x = NULL,y = paste0(Cell_type[i],": log10(RPM) class: ",class[l]))+
        stat_summary(fun="mean",color="black")+
        geom_jitter(shape=16,position = position_jitter(0.1),color="black",size=3)+
        theme_classic2()+
        theme(axis.title.x = element_blank(),
              axis.ticks.x=element_blank(),
              axis.title.y = element_blank(),
              axis.text.x = element_blank(),
              axis.text.y = element_blank(),
              legend.text = element_blank(),
              legend.title = element_blank(),
              legend.position = "none")
      print(pt)
      dev.off()
      # ggsave(paste0("./plot_labs/class/",class[l],"/",geneset_used,"/",as.character(Cell_type[i]),".png"),width = 17.6*0.0394, height = 23.5*0.0394)
    }
  }
  
  ############## summary plot
  my_Plot_class <- function(cell_type){
    tmp<-pseudo_class_metadata[which(pseudo_class_metadata$Cell_type==cell_type),]
    ggplot(tmp,aes(x=Inflammatory_Status,y=log10(tmp[,classes]),color=Inflammatory_Status))+
      stat_boxplot(geom = "errorbar",width=0.2,aes(color=Inflammatory_Status))+
      geom_boxplot(notch = F,width=0.6,outlier.shape = NA)+
      scale_color_manual(values=c("#999999","#1C769F"))+
      labs(x = NULL,
           # title = paste0(cell_type),
           #y =paste0("log10(RPM)")
           y = NULL)+
      stat_summary(fun="mean",color="black")+
      geom_jitter(shape=16, position = position_jitter(0.1),color="black",size=0.5)+
      theme_classic2()+
      theme(axis.title.x = element_blank(),
            axis.ticks.x=element_blank(),
            axis.title.y = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_text(face = "bold", family = "Helvetica", size = "5"),
            legend.text = element_blank(),
            legend.title = element_blank(),
            legend.position = "none")
  }
  
  for (l in 1:length(class)) {
    classes<-class[l]
    pdf(file =paste0("./plot_labs/class/Summarised_",geneset_used,"_",classes,".pdf"),width = 88*0.0394,height=117.3*0.0394)
    my_Plot_List_class <- lapply(Cell_type,my_Plot_class)
    pt <- plot_grid(plotlist = my_Plot_List_class, align = "hv", 
                  nrow = 5,ncol = 5, scale = 0.8)
    print(pt)
    dev.off()
  }
  
  # Single cell type multiple class plot
  for (l in 1:length(class)) {
    classes<-class[l]
    pdf(file =paste0("./plot_labs/class/Summarised_",geneset_used,"_",classes,".pdf"),width = 88*0.0394,height=117.3*0.0394)
    my_Plot_List_class <- lapply(Cell_type,my_Plot_class)
    pt <- plot_grid(plotlist = my_Plot_List_class, align = "hv", 
                    nrow = 5,ncol = 5, scale = 0.8)
    print(pt)
    dev.off()
  }
}

######################################## Supercentenarian cohort ##################################
setwd("../Data/Single_Cell/Supercentenarian")
library("data.table")
library("Matrix")
library("Seurat")

mat <- fread("../Data/Single_Cell/Supercentenarian/01.UMI.txt")
setDF(mat)
# Set the barcodes as the row names.
rownames(mat) <- mat[[1]]
mat[[1]] <- NULL
# Transpose the data and convert to sparse matrix.
mat <- as(t(as.matrix(mat)), "sparseMatrix")
# Create the Seurat object.
seu <- CreateSeuratObject(counts=t(mat))
Cell_Barcodes <- read_delim("Supercentenarian/03.Cell.Barcodes.txt", 
                            delim = "\t", escape_double = FALSE, 
                            col_names = FALSE, trim_ws = TRUE)
colnames(Cell_Barcodes)<-c("cell_barcode","patient","group")
Cell_Barcodes<-as.data.frame(Cell_Barcodes)
rownames(Cell_Barcodes)<-Cell_Barcodes$cell_barcode
Cell_Barcodes<-Cell_Barcodes[colnames(as.matrix(seu@assays$RNA@counts)),]
seu<-AddMetaData(seu,Cell_Barcodes)
seu@meta.data
#saveRDS(seu,"./Supercentenarian/Supercentenarian_scRNA.rds")

seu<-readRDS("Supercentenarian_scRNA.rds")
View(seu@meta.data)
clusters <- read_csv("clusters.csv")
clusters<-as.data.frame(clusters)
rownames(clusters)<-clusters$Barcode
clusters$cell_annotation<-NA
proposed_anno<-c("TC2","NK","TC1","M14","BC","EC", "M16", "MKI", "DC","MGK")
for (i in unique(clusters$Cluster)) {
  clusters[which(clusters$Cluster==i),]$cell_annotation<-proposed_anno[i]
}
seu<-AddMetaData(seu,clusters)
seu@meta.data$cell_barcode
View(seu@meta.data)
#saveRDS(seu,"Supercentenarian_scRNA.rds")

pseudo_class<-readRDS("./pseudo_class.Rds")
pseudo_fam<-readRDS("./pseudo_fam.Rds")
TE_masker <- read.csv("../Data/Single_Cell/hg38_repeat_masker.gz", header=TRUE, sep="\t")
TE_masker <- unique(TE_masker[,c("repName","repClass","repFamily")])
Class_for_comparison<-unique(TE_masker$repClass)
fam_for_comparison<-unique(TE_masker$repFamily)
#View(pseudo_class@meta.data)
#View(pseudo_fam@meta.data)

groupA<-c("SC1","SC2","SC3","SC4","SC5","SC6","SC7")
groupB<-c("CT2", "CT3", "CT4", "CT5")


pseudo_class$group<-NA
pseudo_class@meta.data[which(pseudo_class@meta.data$Sample_ID%in%groupA),]$group<-"Supercentenarian"
pseudo_class@meta.data[which(pseudo_class@meta.data$Sample_ID%in%groupB),]$group<-"Control"
pseudo_class$Status<-NA
pseudo_class@meta.data$Status<-factor(pseudo_class@meta.data$group,levels = c("Control","Supercentenarian"))

pseudo_fam$group<-NA
pseudo_fam@meta.data[which(pseudo_fam@meta.data$Sample_ID%in%groupA),]$group<-"Supercentenarian"
pseudo_fam@meta.data[which(pseudo_fam@meta.data$Sample_ID%in%groupB),]$group<-"Control"
pseudo_fam$Status<-NA
pseudo_fam@meta.data$Status<-factor(pseudo_fam@meta.data$group,levels = c("Control","Supercentenarian"))


pseudo_class_matrix<-as.matrix(pseudo_class@assays$RNA@counts)
pseudo_class_metadata<-pseudo_class@meta.data
pseudo_class_metadata<-pseudo_class_metadata[colnames(pseudo_class_matrix),]
pseudo_class_matrix_normalized<-pseudo_class_matrix[rownames(pseudo_class_matrix)[which(rownames(pseudo_class_matrix)%in%Class_for_comparison)],]
pseudo_class_matrix_normalized_done<-matrix(data = NA,nrow = nrow(pseudo_class_matrix_normalized),ncol = ncol(pseudo_class_matrix_normalized))
rownames(pseudo_class_matrix_normalized_done)<-rownames(pseudo_class_matrix_normalized)
colnames(pseudo_class_matrix_normalized_done)<-colnames(pseudo_class_matrix_normalized)
#delete<-data.frame(sum=colSums(pseudo_class_matrix_normalized))
for (i in 1:ncol(pseudo_class_matrix_normalized)) {
  pseudo_class_matrix_normalized_done[,i]<-
    1000000*pseudo_class_matrix_normalized[,i]/(pseudo_class_metadata[i,]$nCount_RNA)
  #-as.numeric(delete[i,]))
}

pseudo_class_matrix_normalized_done<-t(pseudo_class_matrix_normalized_done)
pseudo_class_metadata<-pseudo_class_metadata[rownames(pseudo_class_matrix_normalized_done),]
pseudo_class_metadata<-cbind(pseudo_class_metadata,pseudo_class_matrix_normalized_done)
pseudo_class_metadata<-pseudo_class_metadata[which(pseudo_class_metadata$Cell_type!="NA"),]
Cell_type<-unique(pseudo_class_metadata$Cell_type)
class<-unique(rownames(pseudo_class_matrix_normalized))

for (i in 1:length(Cell_type)) {
  for (l in 1:length(class)) {
    tmp<-pseudo_class_metadata[which(pseudo_class_metadata$Cell_type==Cell_type[i]),]
    ggplot(tmp,aes(x=Status,y=log10(tmp[,class[l]]),color=Status))+
      stat_boxplot(geom = "errorbar",width=0.2,aes(color=Status))+
      geom_boxplot(notch = F,width=0.6,outlier.shape = NA)+
      scale_color_manual(values=c("#999999","blue"))+
      labs(x = NULL,y = paste0(Cell_type[i],": log10(RPM) class: ",class[l]))+
      #stat_summary(fun="mean",color="black")+
      geom_jitter(shape=16,position = position_jitter(0.1),color="black",size=3)+
      #stat_compare_means()+
      theme_classic2()+
      theme(axis.title.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.y = element_text(size = 16, face = "bold"),
            axis.text.y = element_text(size = 16),
            legend.text = element_blank(),
            legend.title = element_blank())
    ggsave(filename =paste0("./plot/class/",class[l],"/",as.character(Cell_type[i]),".png"),width = 6,height=8)
  }
}

################################### fam ####################################
pseudo_fam_matrix<-as.matrix(pseudo_fam@assays$RNA@counts)
pseudo_fam_metadata<-pseudo_fam@meta.data
pseudo_fam_metadata<-pseudo_fam_metadata[colnames(pseudo_fam_matrix),]
pseudo_fam_matrix_normalized<-pseudo_fam_matrix[rownames(pseudo_fam_matrix)[which(rownames(pseudo_fam_matrix)%in%fam_for_comparison)],]
pseudo_fam_matrix_normalized_done<-matrix(data = NA,nrow = nrow(pseudo_fam_matrix_normalized),ncol = ncol(pseudo_fam_matrix_normalized))
rownames(pseudo_fam_matrix_normalized_done)<-rownames(pseudo_fam_matrix_normalized)
colnames(pseudo_fam_matrix_normalized_done)<-colnames(pseudo_fam_matrix_normalized)
#delete<-data.frame(sum=colSums(pseudo_fam_matrix_normalized))
for (i in 1:ncol(pseudo_fam_matrix_normalized)) {
  pseudo_fam_matrix_normalized_done[,i]<-
    1000000*pseudo_fam_matrix_normalized[,i]/(pseudo_fam_metadata[i,]$nCount_RNA)
  #-as.numeric(delete[i,]))
}
pseudo_fam_matrix_normalized_done<-t(pseudo_fam_matrix_normalized_done)
pseudo_fam_metadata<-pseudo_fam_metadata[rownames(pseudo_fam_matrix_normalized_done),]
pseudo_fam_metadata<-cbind(pseudo_fam_metadata,pseudo_fam_matrix_normalized_done)
pseudo_fam_metadata<-pseudo_fam_metadata[which(pseudo_fam_metadata$Cell_type!="NA"),]
Cell_type<-unique(pseudo_fam_metadata$Cell_type)
fam<-unique(rownames(pseudo_fam_matrix_normalized))
for (i in 1:length(Cell_type)) {
  for (l in 1:length(fam)) {
    tmp<-pseudo_fam_metadata[which(pseudo_fam_metadata$Cell_type==Cell_type[i]),]
    ggplot(tmp,aes(x=Status,y=log10(tmp[,fam[l]]),color=Status))+
      stat_boxplot(geom = "errorbar",width=0.2,aes(color=Status))+
      geom_boxplot(notch = F,width=0.6,outlier.shape = NA)+
      scale_color_manual(values=c("#999999","blue"))+
      labs(x = NULL,y = paste0(Cell_type[i],": log10(RPM) fam: ",fam[l]))+
      #stat_summary(fun="mean",color="black")+
      geom_jitter(shape=16, position = position_jitter(0.1),color="black",size=3)+
      #stat_compare_means()+
      theme_classic2()+
      theme(axis.title.x = element_blank(),
            axis.ticks.x=element_blank(),
            axis.title.y = element_text(size = 16, face = "bold"),
            axis.text.y = element_text(size = 16),
            legend.text = element_blank(),
            legend.title = element_blank())
    ggsave(filename =paste0("./plot/fam/",fam[l],"/",as.character(Cell_type[i]),".png"),width = 6,height=8)
  }
}


################ summary plot
my_Plot_class <- function(cell_type){
  tmp<-pseudo_class_metadata[which(pseudo_class_metadata$Cell_type==cell_type),]
  ggplot(tmp,aes(x=Status,y=log10(tmp[,classes]),color=Status))+
    stat_boxplot(geom = "errorbar",width=0.2,aes(color=Status))+
    geom_boxplot(notch = F,width=0.6,outlier.shape = NA)+
    scale_color_manual(values=c("#999999","blue"))+
    labs(x = NULL,
         title = paste0(cell_type),
         #y =paste0("log10(RPM)")
         y =NULL)+
    stat_summary(fun="mean",color="black")+
    geom_jitter(shape=16, position = position_jitter(0.1),color="black",size=3)+
    stat_compare_means()+
    theme_classic2()+
    theme(axis.title.x = element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.x=element_blank(),
          axis.title.y = element_text(size = 12, face = "bold"),
          axis.text.y = element_text(size = 12),
          legend.text = element_blank(),
          legend.title = element_blank(),
          legend.position = "none",
          plot.title = element_text(size = 14, color = "black",face = "bold"))
}

my_Plot_fam <- function(cell_type){
  tmp<-pseudo_fam_metadata[which(pseudo_fam_metadata$Cell_type==cell_type),]
  ggplot(tmp,aes(x=Status,y=log10(tmp[,family]),color=Status))+
    stat_boxplot(geom = "errorbar",width=0.2,aes(color=Status))+
    geom_boxplot(notch = F,width=0.6,outlier.shape = NA)+
    scale_color_manual(values=c("#999999","blue"))+
    labs(x = NULL,
         title = paste0(cell_type),
         #y =paste0("log10(RPM)")
         y =NULL)+
    stat_summary(fun="mean",color="black")+
    geom_jitter(shape=16, position = position_jitter(0.1),color="black",size=3)+
    stat_compare_means()+
    theme_classic2()+
    theme(axis.title.x = element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.x=element_blank(),
          axis.title.y = element_text(size = 12, face = "bold"),
          axis.text.y = element_text(size = 12),
          legend.text = element_blank(),
          legend.title = element_blank(),
          legend.position = "none",
          plot.title = element_text(size = 14, color = "black",face = "bold"))
}

Cell_type<-unique(pseudo_class_metadata$Cell_type)
for (l in 1:length(class)) {
  classes<-class[l]
  pdf(file =paste0("./plot/class/Summarised_",classes,".pdf"),width = 20,height=15)
  my_Plot_List_class <- lapply(Cell_type,my_Plot_class)
  pt<-plot_grid(plotlist = my_Plot_List_class, align = "hv", 
                nrow = 2,ncol = 5)
  print(pt)
  dev.off()
}

Cell_type<-unique(pseudo_fam_metadata$Cell_type)
for (l in 1:length(fam)) {
  family<-fam[l]
  pdf(file =paste0("./plot/fam/Summarised_",family,".pdf"),width = 20,height=15)
  my_Plot_List_fam <- lapply(Cell_type,my_Plot_fam)
  pt<-plot_grid(plotlist = my_Plot_List_fam, align = "hv", 
                nrow = 2,ncol = 5)
  print(pt)
  dev.off()
}

########################### Inflammation ########################
setwd("../Data/Single_Cell/Supercentenarian/Inflammatory_analysis/")
pseudo_scRNA<-readRDS("../Data/Single_Cell/Supercentenarian/pseudo_class.Rds")
#pseudo_scATAC<-readRDS("../Data/Single_Cell/scATAC-seq_PBMC/pseudo_class.rds")
TE_masker <- read.csv("../Data/Single_Cell/hg38_repeat_masker.gz", header=TRUE, sep="\t")
TE_masker <- unique(TE_masker[,c("repName","repClass","repFamily")])
Class_for_comparison<-unique(TE_masker$repClass)

################## inflammation genelists ######################
Inflamm_gene_list <- read_excel("../Data/Single_Cell/Inflammatory_analysis/Supplementary_tables (1).xlsx",sheet = "S1")
Inflamm_gene_list<-as.data.frame(Inflamm_gene_list)
Lists<-colnames(Inflamm_gene_list)
pseudo_scRNA_count<-as.matrix(pseudo_scRNA@assays$RNA@counts)
library(org.Hs.eg.db)
Inflamm_gene_list[1,]$`IFN-I`<-Inflamm_gene_list[2,]$`IFN-I`
Inflamm_gene_list_2<-Inflamm_gene_list
Inflamm_gene_list_2$`IFN-I` <-mapIds(org.Hs.eg.db, keys = as.character(Inflamm_gene_list$`IFN-I`), column="ENSEMBL", keytype = "SYMBOL")
Inflamm_gene_list_2$`IFN-I` <-Inflamm_gene_list_2$`IFN-I`
Inflamm_gene_list_2$Senescence<-mapIds(org.Hs.eg.db, keys = as.character(Inflamm_gene_list$Senescence), column="ENSEMBL", keytype = "SYMBOL")
Inflamm_gene_list_2$`Inflammatory chemokines`<-mapIds(org.Hs.eg.db, keys = as.character(Inflamm_gene_list$`Inflammatory chemokines`), column="ENSEMBL", keytype = "SYMBOL")
Inflamm_gene_list_2$Inflammaging<-mapIds(org.Hs.eg.db, keys = as.character(Inflamm_gene_list$Inflammaging), column="ENSEMBL", keytype = "SYMBOL")
Inflamm_gene_list_2$`Inflammatory cytokines`<-mapIds(org.Hs.eg.db, keys = as.character(Inflamm_gene_list$`Inflammatory cytokines`), column="ENSEMBL", keytype = "SYMBOL")
Inflamm_gene_list_2$SASP<-mapIds(org.Hs.eg.db, keys = as.character(Inflamm_gene_list$SASP), column="ENSEMBL", keytype = "SYMBOL")
Inflamm_gene_list<-Inflamm_gene_list_2
subset.matrix <- as.matrix(pseudo_scRNA@assays$RNA@counts)[which(rownames(as.matrix(pseudo_scRNA@assays$RNA@counts))%notin%c("LINE","SINE","LTR")), ] 

for (z in 1:length(Lists)) {
  geneset_used<-Lists[z]
  print(Lists[z])
  Inflamm_genes<-unique(as.character(Inflamm_gene_list[,geneset_used]))[which(!is.na(unique(as.character(Inflamm_gene_list[,geneset_used]))))]
  # Pull the raw expression matrix from the original Seurat object containing only the genes of interest
  object2 <- CreateSeuratObject(subset.matrix) # Create a new Seurat object with just the genes of interest
  orig.ident <- pseudo_scRNA@meta.data # Pull the identities from the original Seurat object as a data.frame
  object2 <- AddMetaData(object = object2, metadata = orig.ident) # Add the idents to the meta.data slot
  # Assign identities for the new Seurat object
  Delete<-colSums(as.matrix(pseudo_scRNA@assays$RNA@counts)[c("SINE","LINE","LTR"),])
  object2$nCount_RNA<-object2$nCount_RNA-Delete
  Matrix_patient<-matrix(NA,nrow = nrow(subset.matrix),ncol = length(unique(object2@meta.data$Sample_ID)))
  rownames(Matrix_patient)<-rownames(subset.matrix)
  colnames(Matrix_patient)<-unique(object2@meta.data$Sample_ID)
  Patient_meta<-data.frame(sample=unique(object2@meta.data$Sample_ID),ncount=NA)
  
  for (i in unique(object2@meta.data$Sample_ID)) {
    ident<-rownames(object2@meta.data[which(object2@meta.data$Sample_ID==i),])
    Matrix_patient[,i]<-rowSums(subset.matrix[,ident])
    Patient_meta[which(Patient_meta$sample==i),]$ncount<-sum(subset.matrix[,ident])
  }
  
  Matrix_patient_normalized_done<-matrix(data = NA,nrow = nrow(Matrix_patient),ncol = ncol(Matrix_patient))
  rownames(Matrix_patient_normalized_done)<-rownames(Matrix_patient)
  colnames(Matrix_patient_normalized_done)<-colnames(Matrix_patient)
  
  for (i in 1:ncol(Matrix_patient)) {
    Matrix_patient_normalized_done[,i]<-
      1000000*Matrix_patient[,i]/(Patient_meta[i,]$ncount)
  }
  
  rankData<-rankGenes(Matrix_patient_normalized_done)
  tgfb_gs_up<-as.character(Inflamm_genes)
  scoredf <- simpleScore(rankData, upSet = tgfb_gs_up,knownDirection = T,centerScore = F,dispersionFun = mad)
  High_inflammatory<-rownames(scoredf[which(scoredf$TotalScore>=median(scoredf$TotalScore)),])
  Low_inflammatory<-rownames(scoredf[which(scoredf$TotalScore<median(scoredf$TotalScore)),])
  paste0("High: ",High_inflammatory)
  paste0("Low: ",Low_inflammatory)
  pseudo_class<-readRDS("../Data/Single_Cell/Supercentenarian/pseudo_class.rds")
  pseudo_fam<-readRDS("../Data/Single_Cell/Supercentenarian/pseudo_fam.Rds")
  TE_masker <- read.csv("../Data/Single_Cell/hg38_repeat_masker.gz", header=TRUE, sep="\t")
  TE_masker <- unique(TE_masker[,c("repName","repClass","repFamily")])
  Class_for_comparison<-unique(TE_masker$repClass)
  fam_for_comparison<-unique(TE_masker$repFamily)
  pseudo_class_matrix<-as.matrix(pseudo_class@assays$RNA@counts)
  pseudo_class_metadata<-pseudo_class@meta.data
  pseudo_class_metadata<-pseudo_class_metadata[colnames(pseudo_class_matrix),]
  pseudo_class_matrix_normalized<-pseudo_class_matrix[rownames(pseudo_class_matrix)[which(rownames(pseudo_class_matrix)%in%Class_for_comparison)],]
  pseudo_class_matrix_normalized_done<-matrix(data = NA,nrow = nrow(pseudo_class_matrix_normalized),ncol = ncol(pseudo_class_matrix_normalized))
  rownames(pseudo_class_matrix_normalized_done)<-rownames(pseudo_class_matrix_normalized)
  colnames(pseudo_class_matrix_normalized_done)<-colnames(pseudo_class_matrix_normalized)
  #View(pseudo_class_matrix_normalized)
  #delete<-data.frame(sum=colSums(pseudo_class_matrix_normalized))
  
  for (i in 1:ncol(pseudo_class_matrix_normalized)) {
    pseudo_class_matrix_normalized_done[,i]<-
      1000000*pseudo_class_matrix_normalized[,i]/(pseudo_class_metadata[i,]$nCount_RNA)
    #-as.numeric(delete[i,]))
  }
  pseudo_class_matrix_normalized_done<-t(pseudo_class_matrix_normalized_done)
  pseudo_class_metadata<-pseudo_class_metadata[rownames(pseudo_class_matrix_normalized_done),]
  pseudo_class_metadata<-cbind(pseudo_class_metadata,pseudo_class_matrix_normalized_done)
  pseudo_class_metadata$Status<-factor(pseudo_class_metadata$Sample_ID,levels = c("young","old"))
  groupA<-c("SC1","SC2","SC3","SC4","SC5","SC6","SC7")
  groupB<-c("CT2", "CT3", "CT4", "CT5")
  
  pseudo_class_metadata$Status<-NA
  pseudo_class_metadata[which(pseudo_class_metadata$Sample_ID%in%groupA),]$Status<-"Supercentenarian"
  pseudo_class_metadata[which(pseudo_class_metadata$Sample_ID%in%groupB),]$Status<-"Control"
  pseudo_class_metadata$Status<-factor(pseudo_class_metadata$Status,levels = c("Control","Supercentenarian"))
  pseudo_class_metadata$Inflammatory_Status<-"low_inflammation"
  pseudo_class_metadata[which(pseudo_class_metadata$Sample_ID%in%High_inflammatory),]$Inflammatory_Status<-"high_inflammation"
  pseudo_class_metadata$Inflammatory_Status<-factor(pseudo_class_metadata$Inflammatory_Status,levels = c("low_inflammation","high_inflammation"))
  Cell_type<-unique(pseudo_class_metadata$Cell_type)
  class<-unique(rownames(pseudo_class_matrix_normalized))
  
  for (i in 1:length(Cell_type)) {
    for (l in 1:length(class)) {
      tmp<-pseudo_class_metadata[which(pseudo_class_metadata$Cell_type==Cell_type[i]),]
      ggplot(tmp,aes(x=Inflammatory_Status,y=log10(tmp[,class[l]]),color=Inflammatory_Status))+
        stat_boxplot(geom = "errorbar",width=0.2,aes(color=Inflammatory_Status))+
        geom_boxplot(notch = F,width=0.6,outlier.shape = NA)+
        scale_color_manual(values=c("#999999","#1C769F"))+
        labs(x = NULL,y = paste0(Cell_type[i],": log10(RPM) class: ",class[l]))+
        stat_summary(fun="mean",color="black")+
        geom_jitter(shape=16,position = position_jitter(0.1),color="black",size=3)+
        stat_compare_means()+
        theme_classic2()+
        theme(axis.title.x = element_blank(),
              axis.ticks.x=element_blank(),
              axis.title.y = element_text(size = 16, face = "bold"),
              axis.text.y = element_text(size = 16),
              legend.text = element_blank(),
              legend.title = element_blank())
      ggsave(paste0("./",class[l],"/",geneset_used,"/",as.character(Cell_type[i]),".png"),width = 6,height=8)
    }
  }
  
  ############## summary plot
  my_Plot_class <- function(cell_type){
    tmp<-pseudo_class_metadata[which(pseudo_class_metadata$Cell_type==cell_type),]
    ggplot(tmp,aes(x=Inflammatory_Status,y=log10(tmp[,classes]),color=Inflammatory_Status))+
      stat_boxplot(geom = "errorbar",width=0.2,aes(color=Inflammatory_Status))+
      geom_boxplot(notch = F,width=0.6,outlier.shape = NA)+
      scale_color_manual(values=c("#999999","#1C769F"))+
      labs(x = NULL,
           title = paste0(cell_type),
           #y =paste0("log10(RPM)")
           y =NULL)+
      stat_summary(fun="mean",color="black")+
      geom_jitter(shape=16, position = position_jitter(0.1),color="black",size=3)+
      stat_compare_means()+theme_classic2()+
      theme(axis.title.x = element_blank(),
            axis.ticks.x=element_blank(),
            axis.text.x=element_blank(),
            axis.title.y = element_text(size = 12, face = "bold"),
            axis.text.y = element_text(size = 12),
            legend.text = element_blank(),
            legend.title = element_blank(),
            legend.position = "none",
            plot.title = element_text(size = 14, color = "black",face = "bold"))
  }
  
  for (l in 1:length(class)) {
    classes<-class[l]
    pdf(file =paste0("./Summarised_",geneset_used,"_",classes,".pdf"),width = 15,height=10)
    my_Plot_List_class <- lapply(Cell_type,my_Plot_class)
    pt<-plot_grid(plotlist = my_Plot_List_class, align = "hv", 
                  nrow = 2,ncol = 5)
    print(pt)
    dev.off()
  }
}