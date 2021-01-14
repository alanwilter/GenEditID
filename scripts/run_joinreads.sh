#!/bin/bash
#SBATCH --partition=general
#SBATCH --mem=4096
#SBATCH --job-name joinreads
#SBATCH --output joinreads.out

echo 'Joining paired reads'

echo `pwd`

for f in *.s_1.r_1.fq.gz
  do
    echo $f
    pear -f $f -r `echo $f | cut -d'.' -f1-4`.r_2.fq.gz -o `echo $f | cut -d'.' -f1-4`
  done

for f in *.assembled.fastq
  do
    gzip $f
    mv $f.gz `echo $f | cut -d'.' -f1-4`.fqjoin.gz
  done

rm *.discarded.* *.unassembled.*

echo 'Job done'
