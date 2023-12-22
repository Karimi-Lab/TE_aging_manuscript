# Aging manuscript (to be updated with actual title)
Scripts and data used in the manuscript. 

<img width="1246" alt="Screenshot 2023-12-15 at 15 31 46" src="https://github.com/Karimi-Lab/TE_aging_manuscript/assets/25477052/03ae60ec-0969-4aff-b3be-876b67948210">

## How to use the scripts
_generate_singscores.R_: Steps needed to generate singscsores from Aging-related gene sets using Illumina HumanHT-12 v4. <br>
_methylation_avail_probes.R_: Check Illumina Infinium Human-Methylation probe set in used methylation datasets. Used to generate _Supplementary Figure 1b_. <br>
_scRNA_analysis.R_: Steps needed to run the scRNA analysis and generate related figures. <br>
_plots.Rmd_: Steps needed to run the statisctical analyses and generate related plots and figures. <br><br>

The scripts necessary to convert scRNA data to pseudo-bulk RNA for RTE classes and families can be found under _Scripts/scRNA_Pseudobulk_. 

### Generate gene expression probe availability bars (_Supp. Fig. 1a_)

<img width="628" alt="Screenshot 2023-12-22 at 10 27 02" src="https://github.com/Karimi-Lab/TE_aging_manuscript/assets/25477052/50549c21-b9e1-4710-add3-59d234f518ed">

Run the "Probe bar chart" chunk in _plots.Rmd_. The necessary file _"missing_genes_231221.xlsx"_ is created by running the _generate_singscores.R_ script. Legend and labels were added manually. <br>

### Generate boxplots split by gene set score quartiles for different cohorts and TEs (_Supp. Fig. 6-7_)

<img width="611" alt="Screenshot 2023-12-22 at 10 37 17" src="https://github.com/Karimi-Lab/TE_aging_manuscript/assets/25477052/ecd0bc6d-036c-4090-920f-7d879c3c75a7">

The _"quartile.all"_ function takes cohort and TE class/family as parameters (_plots.Rmd_). For example, to generate the plot for GARP cohort (GSE48556) and TE class LTR, use the following code.

```
quartile.all(cohort = "GSE48556", te_class = "LTR")
```

### Generate boxplots comparing Control vs. Centenarian groups for gene set scores (_Supp. Fig. 13-15_)

<img width="644" alt="Screenshot 2023-12-22 at 12 39 09" src="https://github.com/Karimi-Lab/TE_aging_manuscript/assets/25477052/2245907e-b46e-4601-acc5-8e191ee0fdbd">

The "Supercentenarian cohort" part in _scRNA_analysis.R_. Creates both individual plots for cell types and a "summary" plot combining all cell types in a single pdf (Summarised_<TE>.pdf). Individual plots will be found in Data > Single_Cell > Supercentenarian > Inflammatory_analysis > plot > class/fam > _TE Class_ > _Cell Type_.pdf <br>
If you encounter an error when saving the file, you may need to create the appropriate subfolders manually.

### Generate GSVA heatmap different cohorts and TEs (_Fig. 3a_)

<img width="681" alt="Screenshot 2023-12-22 at 11 31 46" src="https://github.com/Karimi-Lab/TE_aging_manuscript/assets/25477052/03673c5a-d021-4e6c-b331-09285fab4374">

The _"gsva.heatmap"_ function takes a cohort list and TE class/family list as parameters (_plots.Rmd_).

```
c("LINE", "L1", "L2",
  "LTR", "ERV1", "ERVL", "ERVL-MaLR","ERVK",
  "SINE", "Alu", "MIR")

gsva.heatmap(cohorts_list = c("GSE56045", "GSE48556", "GSE58137"), te_list = all_tes)
```

### Generate cell-type specific boxplots of RTE expression vs. age for PBMC scRNA-seq cohorts (_Fig. 5a-f_)

<img width="675" alt="Screenshot 2023-12-22 at 11 01 14" src="https://github.com/Karimi-Lab/TE_aging_manuscript/assets/25477052/9a534d2d-6133-463d-b048-35e819ad0ed8">

The "Inflammation" part in _scRNA_analysis.R_. Creates both individual plots for cell types and a "summary" plot combining all cell types in a single pdf (Summarised_<Gene_Set>_<TE>.pdf). Individual plots will be found in Data > Single_Cell > Inflammatory_analysis > plot > class/fam > _TE Class_ > _Gene Set_ > _Cell Type_.pdf <br>
If you encounter an error when saving the file, you may need to create the appropriate subfolders manually.

## Additional data required (the folder the file should be put into)
MESA cohort gene expression RDS file (Gene_Expression): https://emckclac-my.sharepoint.com/:u:/g/personal/k2140993_kcl_ac_uk/EePfklf6BYNAiPh_xrss5VkBJsdKhhYBU3zHiNpZ-kvYNQ?e=x51Fa4 <br><br>
hg38 RepeatMasker file (Single_Cell): https://emckclac-my.sharepoint.com/:u:/g/personal/k2140993_kcl_ac_uk/EZFfPw8xHllBs3-5flzFExUBpZtOdGs5L_CS959mVZ5aaw?e=SiRRh2 <br><br>
Japanese cohort scRNA expression RDS file (Single_Cell/Supercentenarian): https://emckclac-my.sharepoint.com/:u:/g/personal/k2140993_kcl_ac_uk/EfHOPUPCjdRDtLARUqNcgecBsD5b1UV6BgN_Zp51vbmYww?e=76txII <br><br>
UMI expression matrix (Single_Cell/Supercentenarian): https://emckclac-my.sharepoint.com/:t:/g/personal/k2140993_kcl_ac_uk/EaF0O8LoEcFCtT2TJ3pnv98BYSBBYLN7kGrRkWRAgaUwwQ?e=TqFpRb 
