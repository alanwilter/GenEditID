#!/bin/bash
#SBATCH --partition=general
#SBATCH --mem=4196
#SBATCH --job-name amplicon-vardict
#SBATCH --output amplicon-vardict.out

export JAVA_OPTS="-Xms1g -Xmx4g"
/home/bioinformatics/software/pipelines/ampliconseq/ampliconseq-pipeline-0.4.1/bin/run-pipeline config.vardict.xml
