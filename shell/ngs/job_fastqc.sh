#!/bin/bash
#SBATCH --partition=general
#SBATCH --mem=4096
#SBATCH --job-name fastqc
#SBATCH --output fastqc.out

echo 'Running FastQC'

mkdir fastqc

cd fastq
echo `pwd`

for f in *.fqjoin.gz
  do
    echo $f
    /home/bioinformatics/pipelinesoftware/genomeediting/fastqc-v0.11.8/fastqc --outdir ../fastqc/ --noextract --format fastq $f
  done

echo 'Job done'
