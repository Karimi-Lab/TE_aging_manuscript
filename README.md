# Aging manuscript (to be updated with actual title)
Scripts and data used in the manuscript. 

<img width="1246" alt="Screenshot 2023-12-15 at 15 31 46" src="https://github.com/Karimi-Lab/TE_aging_manuscript/assets/25477052/03ae60ec-0969-4aff-b3be-876b67948210">

## How to use the scripts
_generate_singscores.R_: Steps needed to generate singscsores from Aging-related gene sets using Illumina HumanHT-12 v4. <br>
_methylation_avail_probes.R_: Check Illumina Infinium Human-Methylation probe set in used methylation datasets. Used to generate _Supplementary Figure 1b_. <br>
_scRNA_analysis.R_: Steps needed to run the scRNA analysis and generate related figures. <br>
_plots.Rmd_: Steps needed to run the statisctical analyses and generate related plots and figures.

### Generate gene expression probe availability bars (_Supp. Fig. 1a_)

<img width="628" alt="Screenshot 2023-12-22 at 10 27 02" src="https://github.com/Karimi-Lab/TE_aging_manuscript/assets/25477052/50549c21-b9e1-4710-add3-59d234f518ed">

Run the "Probe bar chart" chunk. The necessary file _"missing_genes_231221.xlsx"_ is created by running the _generate_singscores.R_ script. <br>

```
cohorts_from <- c("GSE48556", "GSE56045", "GSE58137")
cohorts_to <- c("GARP", "MESA","GTP")

class_probes <- read_excel("../Results/missing_genes_231221.xlsx", sheet = "Class Probes")
colnames(class_probes) <- plyr::mapvalues(colnames(class_probes), 
                                          from = cohorts_from, to = cohorts_to)
class_probes <- class_probes %>% 
  dplyr::filter(Class %in% c("LTR", "LINE", "SINE")) %>%
  dplyr::select(c("Class", "Original", "MESA", "GARP", "GTP"))

class_probes_long <- tidyr::pivot_longer(class_probes, cols = c("Original", "MESA", "GARP", "GTP"))
class_probes_long$Class <- factor(class_probes_long$Class, levels=c("LTR", "LINE", "SINE"))
class_probes_long$name <- factor(class_probes_long$name, levels=c("Original", "MESA", "GARP", "GTP"))

family_probes <- read_excel("../Results/missing_genes_231221.xlsx", sheet = "Family Probes")
colnames(family_probes) <- plyr::mapvalues(colnames(family_probes), 
                                          from = cohorts_from, to = cohorts_to)
family_probes <- family_probes %>% 
  dplyr::filter(Family %in% c("L1", "L2", "ERV1", "ERVL", "ERVL-MaLR","ERVK", "Alu", "MIR")) %>%
  dplyr::select(c("Family", "Original", "MESA", "GARP", "GTP"))

family_probes_long <- tidyr::pivot_longer(family_probes, cols = c("Original", "MESA", "GARP", "GTP"))
family_probes_long$Family <- factor(family_probes_long$Family, levels=c("L1", "L2", "ERV1", "ERVL", 
                                                                "ERVL-MaLR","ERVK", "Alu", "MIR"))
colnames(family_probes_long)[1] <- "Class"

# Merge class and family
merged_bar_long <- rbind(class_probes_long, family_probes_long)

pdf(file = paste0("../Figures/images/bar_class.pdf"),
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
  scale_fill_manual(values = c("Original" = "#20854E", "MESA" = "#F0D600", "GARP" = "#D51416", "GTP" = "#F49200"))
print(g)
dev.off()
```

## Additional data required (the folder the file should be put into)
MESA cohort gene expression RDS file (Gene_Expression): https://emckclac-my.sharepoint.com/:u:/g/personal/k2140993_kcl_ac_uk/EePfklf6BYNAiPh_xrss5VkBJsdKhhYBU3zHiNpZ-kvYNQ?e=x51Fa4 <br><br>
hg38 RepeatMasker file (Single_Cell): https://emckclac-my.sharepoint.com/:u:/g/personal/k2140993_kcl_ac_uk/EZFfPw8xHllBs3-5flzFExUBpZtOdGs5L_CS959mVZ5aaw?e=SiRRh2 <br><br>
Japanese cohort scRNA expression RDS file (Single_Cell/Supercentenarian): https://emckclac-my.sharepoint.com/:u:/g/personal/k2140993_kcl_ac_uk/EfHOPUPCjdRDtLARUqNcgecBsD5b1UV6BgN_Zp51vbmYww?e=76txII <br><br>
UMI expression matrix (Single_Cell/Supercentenarian): https://emckclac-my.sharepoint.com/:t:/g/personal/k2140993_kcl_ac_uk/EaF0O8LoEcFCtT2TJ3pnv98BYSBBYLN7kGrRkWRAgaUwwQ?e=TqFpRb 
