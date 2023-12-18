

TE_masker <- read.csv("/scratch/prj/celgene/CML/scRNASeq-PMBC/scTE/TE_read_CREATE/hg38_repeat_masker.gz", header=TRUE, sep="\t")
TE_masker <- unique(TE_masker[,c("repName","repClass","repFamily")])
LINE <- TE_masker[TE_masker$repClass=="LINE",]
SINE <- TE_masker[TE_masker$repClass=="SINE",]
LTR <- TE_masker[TE_masker$repClass=="LTR",]

dir <-"/scratch/prj/celgene/CML/J-DU000905/J-DU000905/JGAS000230/JGAD000324/bam_scTE/TE_read_CREATE/out_csv_files"

TE_class_read <- function(TE_file){
  
  TE_file_csv <- read.csv(file=TE_file, header=TRUE)
  colnames(TE_file_csv) <- gsub("X.","",colnames(TE_file_csv))
  rownames(TE_file_csv) <- TE_file_csv$barcodes
  TE_file_csv <- TE_file_csv[,-1]
  
  LINE_exp <- TE_file_csv[, colnames(TE_file_csv) %in% LINE$repName]
  LINE_exp <- rowSums (LINE_exp, na.rm = FALSE, dims = 1)
  #size_factor <- rowSums (TE_file_csv, na.rm = FALSE, dims = 1)
  LINE_exp <- as.data.frame(LINE_exp)
  colnames(LINE_exp) <- "LINE"
  # LINE_exp$size_factor <- size_factor
  
  LTR_exp <- TE_file_csv[, colnames(TE_file_csv) %in% LTR$repName]
  LTR_exp <- rowSums (LTR_exp, na.rm = FALSE, dims = 1)
  #size_factor <- rowSums (TE_file_csv, na.rm = FALSE, dims = 1)
  LTR_exp <- as.data.frame(LTR_exp)
  colnames(LTR_exp) <- "LTR"
  # LTR_exp$size_factor <- size_factor
  
  SINE_exp <- TE_file_csv[, colnames(TE_file_csv) %in% SINE$repName]
  SINE_exp <- rowSums (SINE_exp, na.rm = FALSE, dims = 1)
  #size_factor <- rowSums (TE_file_csv, na.rm = FALSE, dims = 1)
  SINE_exp <- as.data.frame(SINE_exp)
  colnames(SINE_exp) <- "SINE"
  # SINE_exp$size_factor <- size_factor
  
  TE_boundle <- cbind(LINE_exp, SINE_exp, LTR_exp)
  
  #saveRDS(TE_boundle, file = sub(".out.csv", "_class_reads.rds", TE_file))
  saveRDS(TE_boundle, file = paste0(TE_file,"_class_reads.rds"))
  
  # return(TE_boundle)
  
}

TE_fam_read <- function(TE_file){
  
  # TE_file <- "CITE_MPAL1_T1.csv"
  #  sample_name <- sub("_out.csv","",TE_file) 
  
  TE_file_csv <- read.csv(file=TE_file, header=TRUE)
  colnames(TE_file_csv) <- gsub("X.","",colnames(TE_file_csv))
  rownames(TE_file_csv) <- TE_file_csv$barcodes
  TE_file_csv <- TE_file_csv[,-1]
  
  LINE_exp <- TE_file_csv[, colnames(TE_file_csv) %in% LINE$repName]
  LINE_exp <- as.data.frame(t(LINE_exp))
  rownames(LINE) <- LINE$repName
  LINE <- LINE[rownames(LINE_exp),]
  LINE_exp$repFamily <- LINE$repFamily
  LINE_exp <- aggregate(. ~ repFamily, data = LINE_exp, sum)
  rownames(LINE_exp) <-  LINE_exp$repFamily
  LINE_exp <- LINE_exp[,-1]
  LINE_exp <- t(LINE_exp)
  #size_factor <- rowSums (TE_file_csv, na.rm = FALSE, dims = 1)
  LINE_exp <- as.data.frame(LINE_exp)
  # LINE_exp$size_factor <- size_factor
  
  SINE_exp <- TE_file_csv[, colnames(TE_file_csv) %in% SINE$repName]
  SINE_exp <- as.data.frame(t(SINE_exp))
  rownames(SINE) <- SINE$repName
  SINE <- SINE[rownames(SINE_exp),]
  SINE_exp$repFamily <- SINE$repFamily
  SINE_exp <- aggregate(. ~ repFamily, data = SINE_exp, sum)
  rownames(SINE_exp) <-  SINE_exp$repFamily
  SINE_exp <- SINE_exp[,-1]
  SINE_exp <- t(SINE_exp)
  #size_factor <- rowSums (TE_file_csv, na.rm = FALSE, dims = 1)
  SINE_exp <- as.data.frame(SINE_exp)
  # SINE_exp$size_factor <- size_factor
  
  LTR_exp <- TE_file_csv[, colnames(TE_file_csv) %in% LTR$repName]
  LTR_exp <- as.data.frame(t(LTR_exp))
  LTR <-LTR[LTR$repName!="LTR91",]
  rownames(LTR) <- LTR$repName
  LTR <- LTR[rownames(LTR_exp),]
  LTR_exp$repFamily <- LTR$repFamily
  LTR_exp <- aggregate(. ~ repFamily, data = LTR_exp, sum)
  rownames(LTR_exp) <-  LTR_exp$repFamily
  LTR_exp <- LTR_exp[,-1]
  LTR_exp <- t(LTR_exp)
  #size_factor <- rowSums (TE_file_csv, na.rm = FALSE, dims = 1)
  LTR_exp <- as.data.frame(LTR_exp)
  # LTR_exp$size_factor <- size_factor
  
  TE_boundle <- cbind(LINE_exp, SINE_exp, LTR_exp)
  
 # saveRDS(TE_boundle, file = sub(".out.csv", "_fam_reads.rds", TE_file))
  saveRDS(TE_boundle, file = paste0(TE_file,"_fam_reads.rds"))
}

## *_out.csv created by sort_CB.sh and run_scTE.sh

files <- Sys.glob(paste0(dir,"/*out.csv"))

for (i in files){
  print(i)
  TE_class_read(i)
  TE_fam_read(i)
}



