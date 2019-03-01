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
    fastqc --outdir ../fastqc/ --noextract --format fastq $f
  done

echo 'Job done'
