#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem 32G
#SBATCH --output=jfish_%J_stdout.txt
#SBATCH --error=jfish_%J_stderr.txt
#SBATCH --job-name=jfish
# 

bash jellyfish.sh