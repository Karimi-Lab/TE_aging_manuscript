# Expression of Autonomous Retrotransposons Correlates with Biological Aging
Scripts and data used in the "Expression of Autonomous Retrotransposons Correlates with Biological Aging" paper. 

<img width="1246" alt="Screenshot 2023-12-15 at 15 31 46" src="https://github.com/Karimi-Lab/TE_aging_manuscript/assets/25477052/03ae60ec-0969-4aff-b3be-876b67948210">

## How to use the scripts
generate_singscores.R: Steps needed to generate singscsores from Aging-related gene sets using Illumina HumanHT-12 v4. <br>
methylation_avail_probes.R: Check Illumina Infinium Human-Methylation probe set in used methylation datasets. Used to generate _Supplementary Figure 1b_. <br>
scRNA_analysis.R: Steps needed to run the scRNA analysis and generate related figures.
plots.Rmd: Steps needed to run the statisctical analyses and generate related plots and figures.

## Additional data required (the folder the file should be put into)
MESA cohort gene expression RDS file (Gene_Expression): https://emckclac-my.sharepoint.com/:u:/g/personal/k2140993_kcl_ac_uk/EePfklf6BYNAiPh_xrss5VkBJsdKhhYBU3zHiNpZ-kvYNQ?e=x51Fa4 <br><br>
hg38 RepeatMasker file (Single_Cell): https://emckclac-my.sharepoint.com/:u:/g/personal/k2140993_kcl_ac_uk/EZFfPw8xHllBs3-5flzFExUBpZtOdGs5L_CS959mVZ5aaw?e=SiRRh2 <br><br>
Japanese cohort scRNA expression RDS file (Single_Cell/Supercentenarian): https://emckclac-my.sharepoint.com/:u:/g/personal/k2140993_kcl_ac_uk/EfHOPUPCjdRDtLARUqNcgecBsD5b1UV6BgN_Zp51vbmYww?e=76txII <br><br>
UMI expression matrix (Single_Cell/Supercentenarian): https://emckclac-my.sharepoint.com/:t:/g/personal/k2140993_kcl_ac_uk/EaF0O8LoEcFCtT2TJ3pnv98BYSBBYLN7kGrRkWRAgaUwwQ?e=TqFpRb 
