#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem 24G
#SBATCH --output=NOFL_Aug_%J_stdout.txt
#SBATCH --error=NOFL_Aug_%J_stderr.txt
#SBATCH --job-name=NOFL_Aug
# 

bash NOFL_augustus.sh
