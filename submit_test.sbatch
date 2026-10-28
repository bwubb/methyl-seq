#!/bin/bash
#SBATCH --job-name=methyl-nf
#SBATCH --partition=genoa-std-mem
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=3-00:00:00
#SBATCH --output=logs/nextflow_%j.out
#SBATCH --error=logs/nextflow_%j.err

cd /vast/projects/knathans/testicular-germ-ce/031-Meythl-Seq
export PATH="$HOME/bin:$PATH"
source /etc/profile
module load samtools

nextflow run main.nf \
    -profile betty \
    --input samples.test.tsv \
    --outdir results \
    -resume
