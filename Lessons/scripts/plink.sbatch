#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem 20G
#SBATCH --output=plink_%J_stdout.txt
#SBATCH --error=plink_%J_stderr.txt
#SBATCH --job-name=plink
# 

bash plink.sh
