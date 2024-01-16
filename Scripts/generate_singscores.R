# Script to generate immunology gene scores using gene sets from
# van Galen (2019) Supp. Table 3 or Rutella (2022) Supp. Table 4

##### SETUP #####
library(singscore)
library(vissE)
library(readxl)
library(SummarizedExperiment)
library(edgeR)
library(biomaRt)
library(GSEABase)
library(illuminaHumanv4.db)
library(limma)

# Read Biomart annotations to get the gene lengths necessary for FPKM conversion
# hsapiens_gene_ensembl is GRCh38
ensembl <- biomaRt::useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
gene_annotations <- biomaRt::getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol', 'start_position',
                                                  'end_position'),
                                   mart = ensembl)
gene_annotations$length <- gene_annotations$end_position - gene_annotations$start_position - 1

# ##### To check if probe exists for gene in annotation db
# x <- illuminaHumanv4SYMBOL
# # Get the probe identifiers that are mapped to a gene symbol
# mapped_probes <- mappedkeys(x)
# # Convert to a list
# xx <- as.list(x[mapped_probes])
# if(length(xx) > 0) {
#   # Get the SYMBOL for the first five probes
#   xx[1:5]
#   # Get the first one
#   xx[[1]]
# }

# Load gene sets 
read_excel_allsheets <- function(filename, tibble = FALSE) {
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}
gene_sets <- read_excel_allsheets("../Data/gene_signatures_230607.xlsx")

run_multiScore <- function(x, gs){
  gl <- as.list(gs)
  gsl <- list()
  for (i in 1:length(gl)){
    gsl[names(gl)[i]] <- GSEABase::GeneSet(unique(unlist(gl[i])[!is.na(unlist(gl[i]))]), setName = as.character(names(gl)[i]))
  }
  gs_col <- GSEABase::GeneSetCollection(gsl)
  
  scores <- singscore::multiScore(x, upSetColc = gs_col)
  return(scores)
}

missing_genes <- list()
avail_probes <- read.csv("../Data/Illumina-HumanHT-12-v4-Expression-BeadChip-TE-Probes.csv")

for (rds_file in list.files("../Data/Gene_Expression/",
                            pattern = "*.rds")){
  ##### Preprocessing #####
  # Load gene counts
  gse <- readRDS(paste0("../Data/Gene_Expression/", rds_file))
  if (typeof(gse) == "list"){
    if (length(gse) > 1){
      exprs_df <- as.data.frame(gse[[1]]@assayData$exprs)
    } else {
      exprs_df <- as.data.frame(gse@assayData$exprs)
    }
  } else if (typeof(gse) == "double") {
    exprs_df <- as.data.frame(gse)
  } else if (typeof(gse) == "S4") { # expression set already
    exprs_df <- as.data.frame(exprs(gse))
  }
  
  symbol <- mget(rownames(exprs_df), illuminaHumanv4SYMBOL, ifnotfound = NA)
  exprs_df$gene_symbol <- as.character(unlist(symbol))
  # remove probes not found in database
  exprs_df_filt <- as.data.frame(exprs_df[!is.na(exprs_df$gene_symbol),])
  # take AVERAGE expression of probes for a single gene
  exprs_avg <- as.data.frame(limma::avereps(exprs_df_filt, ID=exprs_df_filt$gene_symbol))
  genes <- exprs_avg[,'gene_symbol']
  exprs_avg$gene_symbol <- NULL
  exprs_avg <- as.data.frame(sapply(exprs_avg, as.numeric))
  rownames(exprs_avg) <- genes
  
  # Save expression matrix
  # saveRDS(exprs_avg, file = paste0("./Data/Gene_Expression/", strsplit(rds_file, split = "\\.")[[1]][1], "_avg.Rds"))

  ##### Calculate Immune Scores #####
  ranked_df <- singscore::rankGenes(exprs_avg)
  rownames(ranked_df) <- genes

  # Run singscore for all gene sets for TCGA
  sign_scores <- lapply(names(gene_sets), FUN = function(x) run_multiScore(ranked_df, gene_sets[[x]]))
  names(sign_scores) <- names(gene_sets)

  # Add missing genes in each cohort to missing_genes list
  missing_genes[[strsplit(rds_file, split = "\\.")[[1]][1]]] <- unlist(lapply(names(gene_sets),
                                     FUN = function(x){paste0(x, ":", setdiff(gene_sets[[x]][[1]], rownames(ranked_df)))}))

  # Add available TE probes for each cohort
  avail_probes[,strsplit(rds_file, split = "\\.")[[1]][1]] <- unlist(lapply(avail_probes$probe_ID,
                                                                     FUN = function(x){ifelse(x %in% rownames(exprs_df), 1,0)}))

  ##### Write signiture scores to excel file #####
  # open new excel workbook for writing tables
  OUTPUT <- openxlsx::createWorkbook()

  for (n in names(sign_scores)){
    openxlsx::addWorksheet(OUTPUT, n)
    openxlsx::writeData(OUTPUT, sheet = n, x = sign_scores[[n]]$Scores, rowNames = T)
  }

  openxlsx::saveWorkbook(OUTPUT,
                         paste0("../Results/", strsplit(rds_file, split = "\\.")[[1]][1], "_", format(Sys.Date(), "%y%m%d"), ".xlsx"),
                         overwrite = T)
}

avail_probes <- avail_probes[!duplicated(avail_probes[c("probe_ID","Class","Family")]),]

pr_class <- as.data.frame(lapply(colnames(avail_probes)[7:18],
                                 FUN = function(x){aggregate(avail_probes[,x],
                                                             by = list(Class = avail_probes$Class), FUN = sum)}))
pr_class <- pr_class[,-(grep("Class",colnames(pr_class))[-1])]
colnames(pr_class) <- c("Class", colnames(avail_probes)[7:18])

pr_class$Original <- as.numeric(table(avail_probes$Class))

pr_family <- as.data.frame(lapply(colnames(avail_probes)[7:18],
                                  FUN = function(x){aggregate(avail_probes[,x],
                                                              by = list(Family = avail_probes$Family), FUN = sum)}))
pr_family <- pr_family[,-(grep("Family",colnames(pr_family))[-1])]
colnames(pr_family) <- c("Family", colnames(avail_probes)[7:18])

pr_family$Original <- as.numeric(table(avail_probes$Family) )

OUTPUT <- openxlsx::createWorkbook()

for (n in names(missing_genes)){
  openxlsx::addWorksheet(OUTPUT, n)
  openxlsx::writeData(OUTPUT, sheet = n, x = missing_genes[[n]], rowNames = T)
}

openxlsx::addWorksheet(OUTPUT, "Class Probes")
openxlsx::writeData(OUTPUT, sheet = "Class Probes", x = pr_class, rowNames = F)
openxlsx::addWorksheet(OUTPUT, "Family Probes")
openxlsx::writeData(OUTPUT, sheet = "Family Probes", x = pr_family, rowNames = F)

openxlsx::saveWorkbook(OUTPUT, paste0("../Results/missing_genes_", format(Sys.Date(), "%y%m%d"),".xlsx"), overwrite = T)
