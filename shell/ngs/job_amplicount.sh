#!/bin/bash
#SBATCH --partition=general
#SBATCH --mem=4096
#SBATCH --job-name amplicount
#SBATCH --output amplicount.out

source ~/GenEditID/venv/bin/activate
python ~/GenEditID/python/scripts/run_ampli_count.py --fastqdir=fastq/
