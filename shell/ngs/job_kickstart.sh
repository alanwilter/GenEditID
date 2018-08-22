#!/bin/bash
#SBATCH --partition=general
#SBATCH --mem=4096
#SBATCH --job-name kickstart
#SBATCH --output kickstart.out

/home/bioinformatics/software/pipelines/kickstart/current/bin/kickstart \
  --aligner=bwamem \
  --genome=GRCh38_hs38d1 \
  --library=$1 \
  --fastq-only --no-mark-duplicates
