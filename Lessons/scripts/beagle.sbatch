#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem 20G
#SBATCH --output=beagle_%J_stdout.txt
#SBATCH --error=beagle_%J_stderr.txt
#SBATCH --job-name=beagle
# 

bash beagle.sh
