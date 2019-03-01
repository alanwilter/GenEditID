#!/bin/bash
#SBATCH --partition=general
#SBATCH --mem=4096
#SBATCH --job-name faToTwoBit
#SBATCH --output faToTwoBit.out

faToTwoBit \
  /scratcha/bioinformatics/reference_data/reference_genomes/homo_sapiens/GRCh38_hs38d1/fasta/hsa.GRCh38_hs38d1.fa \
  GRCh38_hs38d1.2bit

#faToTwoBit \
#  /scratcha/bioinformatics/reference_data/reference_genomes/homo_sapiens/GRCh38/fasta/hsa.GRCh38.fa \
#  GRCh38.2bit 
