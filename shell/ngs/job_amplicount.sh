#!/bin/bash
#SBATCH --partition=general
#SBATCH --mem=4096
#SBATCH --job-name amplicount
#SBATCH --output amplicount.out

source /home/pajon01/genome-editing/venv/bin/activate
python /home/pajon01/genome-editing/python/scripts/amplicount.py --fastqdir=fastq/
