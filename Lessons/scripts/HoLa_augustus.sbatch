#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem 24G
#SBATCH --output=HoLa_Aug_%J_stdout.txt
#SBATCH --error=HoLa_Aug_%J_stderr.txt
#SBATCH --job-name=HoLa_Aug
# 

bash HoLa_augustus.sh
