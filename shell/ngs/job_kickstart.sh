#!/bin/bash
#SBATCH --ntasks=1 # Number of cores
#SBATCH --nodes=1 # Ensure that all cores are on one machine
#SBATCH --time=0-00:00 # 1-00:00 Runtime in D-HH:MM
#SBATCH --mem=4000 # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --job-name kickstart-fetch-fastq
#SBATCH --output kickstart-fetch-fastq.%j.out

/home/bioinformatics/software/pipelines/kickstart/current/bin/kickstart \
  --aligner=bwamem \
  --genome=GRCh38 \
  --library=SLX-13775 \
  --fastq-only
