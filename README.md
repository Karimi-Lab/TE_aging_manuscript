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

Run the "Probe bar chart" chunk in _plots.Rmd_. The necessary file _"missing_genes_231221.xlsx"_ is created by running the _generate_singscores.R_ script. Legend and labels were added manually. <br>

### Generate boxplots split by gene set score quartiles for different cohorts and TEs (_Supp. Fig. 6-7_)

<img width="611" alt="Screenshot 2023-12-22 at 10 37 17" src="https://github.com/Karimi-Lab/TE_aging_manuscript/assets/25477052/ecd0bc6d-036c-4090-920f-7d879c3c75a7">

The _"quartile.all"_ function takes cohort and TE class/family as parameters (_plots.Rmd_).

```
quartile.all(cohort = "GSE48556", te_class = "LTR")
```

## Additional data required (the folder the file should be put into)
MESA cohort gene expression RDS file (Gene_Expression): https://emckclac-my.sharepoint.com/:u:/g/personal/k2140993_kcl_ac_uk/EePfklf6BYNAiPh_xrss5VkBJsdKhhYBU3zHiNpZ-kvYNQ?e=x51Fa4 <br><br>
hg38 RepeatMasker file (Single_Cell): https://emckclac-my.sharepoint.com/:u:/g/personal/k2140993_kcl_ac_uk/EZFfPw8xHllBs3-5flzFExUBpZtOdGs5L_CS959mVZ5aaw?e=SiRRh2 <br><br>
Japanese cohort scRNA expression RDS file (Single_Cell/Supercentenarian): https://emckclac-my.sharepoint.com/:u:/g/personal/k2140993_kcl_ac_uk/EfHOPUPCjdRDtLARUqNcgecBsD5b1UV6BgN_Zp51vbmYww?e=76txII <br><br>
UMI expression matrix (Single_Cell/Supercentenarian): https://emckclac-my.sharepoint.com/:t:/g/personal/k2140993_kcl_ac_uk/EaF0O8LoEcFCtT2TJ3pnv98BYSBBYLN7kGrRkWRAgaUwwQ?e=TqFpRb 
