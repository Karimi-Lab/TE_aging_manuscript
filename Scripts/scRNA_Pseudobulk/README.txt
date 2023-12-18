## Alignment and Pseudobulking

1. Raw fastq files for the 10X were aligned using cellranger:

> cellranger_count.sh


2. scTE run over the samples using the run_scTE.sh script:

> run_scTE.sh

N.B. This script runs scTE in a conda environment, which is easily acheived using the included yml file:

> conda env create -f scTE_env.yml


3. TE seurat objects are created from the scTE output using 

> TE_read_CREATE.sh

  ... a wrapper script for TE_read_CREATE.R

N.B. This script, and subsequent scripts, rely on running R with Seurat installed in a conda environment, which is easily acheived using the included yml file:
 
> conda env create -f seurat-env.yml


4. Samples TE Seurat objects for samples into gene Seurat objects

> AddTE.sh

  ... a wrapper script for AddTE.R


5. Individual sample TE Seurat objects are merged into a single Seurat Object

> MergeTE.sh

  .. wrapper script for MergeTE.R


6. single cell data was the pseudobulked on the basis of the cell populations

> PseudoBulk.sh

  ... a wrapper script for Seurat_PseudoBulk.r




