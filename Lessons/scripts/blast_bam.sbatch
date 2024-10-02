#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem 16G
#SBATCH --output=bam_blast_%J_stdout.txt
#SBATCH --error=bam_blast_%J_stderr.txt
#SBATCH --job-name=bam_blast
# 

bash blast_bam.sh
