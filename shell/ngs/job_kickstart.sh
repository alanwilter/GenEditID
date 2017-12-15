#!/bin/bash
#SBATCH --partition=general
#SBATCH --mem=4096
#SBATCH --job-name kickstart
#SBATCH --output kickstart.out

/home/bioinformatics/software/pipelines/kickstart/current/bin/kickstart \
  --aligner=bwamem \
  --species=homo_sapiens \
  --genome-version=GRCh38_hs38d1 \
  --library=$1 \
  --fastq-only
