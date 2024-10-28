#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem 20G
#SBATCH --output=fst_%J_stdout.txt
#SBATCH --error=fst_%J_stderr.txt
#SBATCH --job-name=fst
# 

bash fst.sh
