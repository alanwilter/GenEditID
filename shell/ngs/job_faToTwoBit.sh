#!/bin/bash
#SBATCH --partition=general
#SBATCH --mem=4096
#SBATCH --job-name faToTwoBit
#SBATCH --output faToTwoBit.out

faToTwoBit hsa.GRCh38_hs38d1.fa GRCh38_hs38d1.2bit

echo 'Job done'
