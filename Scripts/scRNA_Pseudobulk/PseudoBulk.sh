#! /usr/bin/bash

#SBATCH --job-name=PseudoBulk
#SBATCH --partition=celgene_cpu 
#SBATCH --time=1-00:00:00 
#SBATCH --nodes=1 
#SBATCH --cpus-per-task=8 
#SBATCH --mem=64G
#SBATCH --chdir=/users/k21216350/CML/scRNASeq-PMBC/scTE
#SBATCH --output=%x_%j.log

# cd /users/k21216350/CML/scRNASeq-PMBC/scTE

echo $SLURM_SUBMIT_DIR

module load anaconda3
CONDA_BASE=$(conda info --base) 
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate seurat-env

Rscript -e 'source("Seurat_PseudoBulk.r")'


