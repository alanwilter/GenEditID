#!/bin/bash
#SBATCH --partition=general
#SBATCH --mem=4096
#SBATCH --job-name amplicount
#SBATCH --output amplicount.out

source ~/editID/venv/bin/activate
python ~/editID/python/scripts/run_ampli_count.py --fastqdir=fastq/
