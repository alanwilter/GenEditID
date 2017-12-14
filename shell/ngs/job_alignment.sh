#!/bin/bash
#SBATCH --partition=general
#SBATCH --mem=4096
#SBATCH --job-name alignment
#SBATCH --output alignment.out

/home/bioinformatics/software/pipelines/alignment/current/bin/run-pipeline --mode=slurm realignment-meta.xml
