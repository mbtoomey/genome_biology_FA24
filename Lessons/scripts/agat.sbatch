#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem 16G
#SBATCH --output=agat_%J_stdout.txt
#SBATCH --error=agat_%J_stderr.txt
#SBATCH --job-name=agat
# 

bash agat.sh