#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem 24G
#SBATCH --output=gemma_%J_stdout.txt
#SBATCH --error=gemma_%J_stderr.txt
#SBATCH --job-name=gemma
# 

bash gemma_lm.sh
