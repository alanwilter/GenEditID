#!/bin/bash
#SBATCH --partition=general
#SBATCH --mem=4096
#SBATCH --job-name multiqc
#SBATCH --output multiqc.out

echo 'Running MultiQC'

mkdir multiqc
source ~/venv-multiqc/bin/activate
multiqc --verbose --dirs --export --outdir multiqc/ fastqc/*fqjoin_fastqc.zip

echo 'Job done'
