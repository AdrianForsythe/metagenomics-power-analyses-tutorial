#!/bin/bash
#SBATCH -A snic2021-5-477
#SBATCH -p node
#SBATCH -N 1
#SBATCH -t 12:00:00
#SBATCH -J power
#SBATCH -o /proj/sllstore2017021/nobackup/ADRIAN/logs/slurm-%A_%a.out
#SBATCH --mail-type=FAIL

source /sw/apps/conda/latest/rackham/etc/profile.d/conda.sh
conda activate /home/adrianf/project-folder/nobackup/ADRIAN/metagenomics-power-analyses-tutorial/r-env

module load R/4.1.1

Rscript /home/adrianf/project-folder/nobackup/ADRIAN/scripts/general-scripts/power-analyses/trial-run-currbio-bear.R

Rscript /home/adrianf/project-folder/nobackup/ADRIAN/scripts/general-scripts/power-analyses/trial-run-full-bear.R
