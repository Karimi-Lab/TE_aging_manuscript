#! /usr/bin/bash

#SBATCH --job-name=Add_Merge_TE
#SBATCH --partition=cpu 
#SBATCH --time=1-00:00:00 
#SBATCH --nodes=1 
#SBATCH --cpus-per-task=8 
#SBATCH --mem=64G 
#SBATCH --output=/scratch/prj/celgene/CML/J-DU000905/J-DU000905/JGAS000230/JGAD000324/bam_scTE/AddTE/Add_TE_%j.log


module load anaconda3
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate seurat-env

cd /scratch/prj/celgene/CML/J-DU000905/J-DU000905/JGAS000230/JGAD000324/bam_scTE/AddTE/

Rscript AddTE.R


