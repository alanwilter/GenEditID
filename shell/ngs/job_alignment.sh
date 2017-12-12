#!/bin/bash
#SBATCH --ntasks=1 # Number of cores
#SBATCH --nodes=1 # Ensure that all cores are on one machine
#SBATCH --time=0-00:00 # 1-00:00 Runtime in D-HH:MM
#SBATCH --mem=4000 # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --job-name alignment-pipeline
#SBATCH --output alignment-pipeline.%j.out

export SOFTWARE_ROOT=/home/bioinformatics/pipelinesoftware/alignment/el7
export JAVA_HOME=/home/bioinformatics/software/java/latest
export JAVA_OPTS=-Xmx1500m
/home/bioinformatics/software/pipelines/alignment/current/bin/run-pipeline --mode=slurm alignment.xml
