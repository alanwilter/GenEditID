#!/bin/bash
#SBATCH --partition=general
#SBATCH --mem=4096
#SBATCH --job-name mergereads
#SBATCH --output mergereads.out

echo 'Merging paired reads'

echo `pwd`

# reverse complement R2
for f in *.s_1.r_2.fq.gz
  do
    echo $f
    # reverse complement R2
    seqkit seq -t dna -v -r -p $f -o `echo $f | cut -d'.' -f1-5`.fqrc.gz
    # combine R1 & R2
    cat `echo $f | cut -d'.' -f1-4`.r_1.fq.gz `echo $f | cut -d'.' -f1-4`.r_2.fqrc.gz | zcat | gzip > `echo $f | cut -d'.' -f1-4`.fqjoin.gz
    rm `echo $f | cut -d'.' -f1-4`.r_2.fqrc.gz
  done

echo 'Job done'
