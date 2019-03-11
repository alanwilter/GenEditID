#!/bin/bash
#SBATCH --partition=general
#SBATCH --mem=4096
#SBATCH --job-name joinreads
#SBATCH --output joinreads.out

echo 'Joining paired reads'

cd fastq
echo `pwd`

for f in *.s_1.r_1.fq.gz
  do
    echo $f
    /home/bioinformatics/pipelinesoftware/genomeediting/fastq-join $f `echo $f | cut -d'.' -f1-4`.r_2.fq.gz -o `echo $f | cut -d'.' -f1-4`.fq
  done

for f in *.fqjoin
  do
    gzip $f
  done

rm *.fqun?

echo 'Job done'
