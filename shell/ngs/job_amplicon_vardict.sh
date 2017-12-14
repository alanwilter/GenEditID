#!/bin/bash
#SBATCH --partition=general
#SBATCH --mem=4096
#SBATCH --job-name amplicon-vardict
#SBATCH --output amplicon-vardict.out

/home/bioinformatics/software/pipelines/ampliconseq/ampliconseq-pipeline-0.4.1/bin/run-pipeline --mode slurm config.vardict.xml
