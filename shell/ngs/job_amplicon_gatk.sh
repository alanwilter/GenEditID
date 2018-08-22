#!/bin/bash
#SBATCH --partition=general
#SBATCH --mem=4196
#SBATCH --job-name amplicon-gatk
#SBATCH --output amplicon-gatk.out

export JAVA_OPTS="-Xms1g -Xmx4g"
/home/bioinformatics/software/pipelines/ampliconseq/ampliconseq-pipeline-0.4.1/bin/run-pipeline config.gatk.xml
