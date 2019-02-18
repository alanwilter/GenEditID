#!/bin/bash
#SBATCH --partition=general
#SBATCH --mem=4096
#SBATCH --job-name readcounts
#SBATCH --output readcounts.out

source /home/pajon01/genome-editing/venv/bin/activate
python /home/pajon01/genome-editing/python/scripts/count_amplicons.py GEP00005_IRX5.yml > GEP00005_counts_IRX5_join.csv
