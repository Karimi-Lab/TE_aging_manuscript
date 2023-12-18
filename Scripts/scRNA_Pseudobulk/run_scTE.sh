#!/bin/bash
#SBATCH --job-name=scTE
#SBATCH --nodes=1
#SBATCH --partition=celgene_cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=2-00:00:00
#SBATCH --output=/scratch/prj/celgene/CML/J-DU000905/J-DU000905/JGAS000230/JGAD000324/bam_scTE/%x_%A-%a.log
#SBATCH --array=12

module load anaconda3
module load samtools
CONDA_BASE=$(conda info --base) 
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate scTE_env
export PATH=$PATH:/scratch/prj/celgene/shared/groups/karimi/Ian_R_Thompson/scTE/bin/

cd /scratch/prj/celgene/CML/J-DU000905/J-DU000905/JGAS000230/JGAD000324/bam_scTE

file=$(ls ./*/*.bam | sed -n ${SLURM_ARRAY_TASK_ID}p | rev | cut -d "/" -f 1 |rev )

samtools view -@ $SLURM_CPUS_PER_TASK -h -d CB -b -o $file ./*/${file}

scTE -i ${file} \
  -p $SLURM_CPUS_PER_TASK \
  -o ${file%.bam}.out \
  --keeptmp True \
  -CB CB -UMI UB \
  -x /scratch/prj/celgene/shared/groups/karimi/Ian_R_Thompson/scTE/hg38.exclusive.idx 


