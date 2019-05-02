#!/bin/bash
#SBATCH --partition=general
#SBATCH --mem=4096
#SBATCH --job-name amplicount
#SBATCH --output amplicount.out

source /home/pajon01/editid/venv/bin/activate
python /home/pajon01/editid/python/scripts/run_ampli_count.py --fastqdir=fastq/
