#! /bin/bash

#SBATCH --output=/scratch/prj/celgene/CML/HRA000203/%x_%A-%a.log
#SBATCH --partition=celgene_cpu
#SBATCH --time=2-00:00:00
#SBATCH --nodes=1
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --job-name=cellranger_scRNA
#SBATCH --array=1-52

cellranger="/scratch/prj/celgene/shared/groups/karimi/Ian_R_Thompson/cellranger-7.1.0/cellranger"
wd="/scratch/prj/celgene/CML/HRA000203/"
ref="/scratch/prj/celgene/shared/groups/karimi/Ian_R_Thompson/refdata-gex-GRCh38-2020-A"
fastq_dir="${wd}/Fastq"

cd $wd

id=$(find ${fastq_dir}/* -type d | sed -n ${SLURM_ARRAY_TASK_ID}p | rev | cut -d "/" -f 1 | rev )

sample=$(find ${fastq_dir}/* -type d | sed -n ${SLURM_ARRAY_TASK_ID}p | rev | cut -d "/" -f 1 | rev )

mem=$(($SLURM_MEM_PER_NODE / 1024))

$cellranger count \
  --id $id \
  --sample $sample \
  --transcriptome $ref \
  --fastqs ${fastq_dir} \
  --localmem $mem  \
  --localcores $SLURM_CPUS_PER_TASK
 
