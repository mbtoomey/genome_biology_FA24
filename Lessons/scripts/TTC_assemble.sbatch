#!/bin/bash
#
#SBATCH --partition=64gb_24core
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem 60G
#SBATCH --output=TTC_%J_stdout.txt
#SBATCH --error=TTC_%J_stderr.txt
#SBATCH --job-name=TTC
# 

bash TTC_assemble.sh