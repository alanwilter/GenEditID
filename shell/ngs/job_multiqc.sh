#!/bin/bash
#SBATCH --partition=general
#SBATCH --mem=4096
#SBATCH --job-name multiqc
#SBATCH --output multiqc.out

echo 'Running MultiQC'

mkdir multiqc
source /home/bioinformatics/pipelinesoftware/genomeediting/venv-multiqc/bin/activate
multiqc --verbose --dirs --export --outdir multiqc/ fastqc/*.fqjoin_fastqc.zip bam/*.bwamem.homo_sapiens.alignment_metrics.txt

echo 'Job done'
