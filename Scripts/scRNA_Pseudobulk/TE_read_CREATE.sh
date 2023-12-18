#! /usr/bin/bash

#SBATCH --job-name=TE_reads
#SBATCH --partition=celgene_cpu 
#SBATCH --time=1-00:00:00 
#SBATCH --nodes=1 
#SBATCH --cpus-per-task=8 
#SBATCH --mem=64G 
#SBATCH --output=/scratch/prj/celgene/CML/J-DU000905/J-DU000905/JGAS000230/JGAD000324/bam_scTE/TE_read_CREATE/TE_read_CREATE_%j.log


module load anaconda3
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate seurat-env

cd /scratch/prj/celgene/CML/J-DU000905/J-DU000905/JGAS000230/JGAD000324/bam_scTE/TE_read_CREATE/

Rscript TE_read_CREATE.R


