#!/bin/bash
#SBATCH --partition general
#SBATCH --ntasks=1 # Number of cores
#SBATCH --nodes=1 # Ensure that all cores are on one machine
#SBATCH --time=0-00:00 # 1-00:00 Runtime in D-HH:MM
#SBATCH --mem=4000 # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --job-name amplicon-pipeline
#SBATCH --output amplicon-pipeline.%j.out

export JAVA_OPTS=-Xmx10240m
perlbrew use perl-5.16.0 && /home/bioinformatics/software/pipelines/ampliconseq/ampliconseq-pipeline-0.4.1/bin/run-pipeline \
  --mode slurm \
  config.gatk.xml
