#!/bin/bash
#SBATCH --partition=general
#SBATCH --mem=4096
#SBATCH --job-name amplicon
#SBATCH --output amplicon.out

/home/bioinformatics/software/pipelines/ampliconseq/ampliconseq-pipeline-0.4.1/bin/run-pipeline --mode slurm config.gatk.xml
