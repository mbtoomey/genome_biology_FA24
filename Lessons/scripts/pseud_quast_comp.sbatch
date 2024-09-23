#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem 16G
#SBATCH --output=quast_%J_stdout.txt
#SBATCH --error=quast_%J_stderr.txt
#SBATCH --job-name=quast
# 

bash pseud_quast_comp.sh