library(readxl)

avail_probes <- read.csv("../Data/Methylation/Infinium-HumanMethylation450-BeadChip-array-Probes.csv")
avail_probes <- avail_probes[!duplicated(avail_probes$TE_Identity),]

for (methyl.cohort in c("GSE56046", "GSE40279", "E-MTAB-7309", "GSE56105")){
  gse <- readRDS(paste0("../Data/Methylation/", methyl.cohort, "_Methylation_TE_Final.rds"))
  exprs_df <- as.data.frame(gse)
  
  # Remove NA rows
  exprs_df <- na.omit(exprs_df)
  
  avail_probes[,methyl.cohort] <- unlist(lapply(avail_probes$TE_Identity,
                                                FUN = function(x){ifelse(x %in% rownames(exprs_df), 1,0)}))
  
}

pr_class <- as.data.frame(lapply(colnames(avail_probes)[7:10],
                                 FUN = function(x){aggregate(avail_probes[,x],
                                                             by = list(Class = avail_probes$Class), FUN = sum)}))
pr_class <- pr_class[,-(grep("Class",colnames(pr_class))[-1])]
colnames(pr_class) <- c("Class", colnames(avail_probes)[7:10])

pr_class$Available <- as.numeric(table(avail_probes$Class))

pr_family <- as.data.frame(lapply(colnames(avail_probes)[7:10],
                                  FUN = function(x){aggregate(avail_probes[,x],
                                                              by = list(Family = avail_probes$Family), FUN = sum)}))
pr_family <- pr_family[,-(grep("Family",colnames(pr_family))[-1])]
colnames(pr_family) <- c("Family", colnames(avail_probes)[7:10])

pr_family$Available <- as.numeric(table(avail_probes$Family))

# Remove ? rows
pr_class <- pr_class[!grepl("\\?",pr_class$Class),]
pr_family <- pr_family[!grepl("\\?",pr_family$Family),]

OUTPUT <- openxlsx::createWorkbook()

openxlsx::addWorksheet(OUTPUT, "Class Probes")
openxlsx::writeData(OUTPUT, sheet = "Class Probes", x = pr_class, rowNames = F)
openxlsx::addWorksheet(OUTPUT, "Family Probes")
openxlsx::writeData(OUTPUT, sheet = "Family Probes", x = pr_family, rowNames = F)

openxlsx::saveWorkbook(OUTPUT, paste0("../Results/methylation_avail_probes_", format(Sys.Date(), "%y%m%d"),".xlsx"), 
                       overwrite = T)

##### Plot #####
class_probes <- read_excel(paste0("../Results/methylation_avail_probes_", format(Sys.Date(), "%y%m%d"),".xlsx"), 
                           sheet = "Class Probes")

class_probes <- class_probes %>% 
  dplyr::filter(Class %in% c("LTR", "LINE", "SINE")) %>%
  dplyr::select(c("Class", "Available", "GSE56046", "GSE40279", "E-MTAB-7309", "GSE56105"))

class_probes_long <- tidyr::pivot_longer(class_probes, cols = c("Available", "GSE56046", "GSE40279", "E-MTAB-7309", "GSE56105"))
class_probes_long$Class <- factor(class_probes_long$Class, levels=c("LTR", "LINE", "SINE"))
class_probes_long$name <- factor(class_probes_long$name, levels=c("Available", "GSE56046", "GSE40279", "E-MTAB-7309", "GSE56105"))

family_probes <- read_excel(paste0("../Results/methylation_avail_probes_", format(Sys.Date(), "%y%m%d"),".xlsx"), 
                            sheet = "Family Probes")

family_probes <- family_probes %>% 
  dplyr::filter(Family %in% c("L1", "L2", "ERV1", "ERVL", "ERVL-MaLR","ERVK", "Alu", "MIR")) %>%
  dplyr::select(c("Family", "Available", "GSE56046", "GSE40279", "E-MTAB-7309", "GSE56105"))

family_probes_long <- tidyr::pivot_longer(family_probes, cols = c("Available", "GSE56046", "GSE40279", "E-MTAB-7309", "GSE56105"))
family_probes_long$Family <- factor(family_probes_long$Family, levels=c("L1", "L2", "ERV1", "ERVL", 
                                                                        "ERVL-MaLR","ERVK", "Alu", "MIR"))
colnames(family_probes_long)[1] <- "Class"

# Merge class and family
merged_bar_long <- rbind(class_probes_long, family_probes_long)

pdf(file = paste0("../Figures/images/methyl_bar_avail.pdf"),
    width = 160*0.0394,
    height = 40*0.0394)
g <- ggplot(data = merged_bar_long, aes(x = Class, y = value, fill = name)) +
  geom_bar(position = "dodge" , stat = "identity") +
  theme_pubr() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 7, family = "Helvetica"),
        legend.position = "none") +
  scale_fill_manual(values = c("Available" = "#20854E", "GSE56046" = "#00A087", "GSE40279" = "#4CBBD5", 
                               "E-MTAB-7309" = "#3C5488", "GSE56105" = "#E64B35"))
print(g)
dev.off()
