#!/bin/bash

if [ $# -lt 4 ]
then
	echo "Usage: $0 <id> <fastq_dir> <queue> <config_file>"
	exit 1
fi

id=$1
fastqdir=$2
queue=$3
configfile=$4

fastqdir=`echo $fastqdir | sed "s|/$||"`

/home/bioinformatics/software/pipelines/alignment/alignment-2.2.2/bin/create-config-file \
	--id=${id} \
	--dataset-id-regex="${id}\.(FLD.*)\..*\.s_[1-8]\.r_[12]\.fq\.gz" \
	--paired-end \
	--platform=ILLUMINA \
	--aligner bwamem \
	--bwa=/home/bioinformatics/software/bwa/bwa-0.7.15/bwa \
	--reference-genome-sequence=/mnt/scratcha/bioinformatics/chilam01/reference_data/homo_sapiens/GRCh38/fasta/hsa.GRCh38.fa \
	--reference-genome-prefix=/mnt/scratcha/bioinformatics/chilam01/reference_data/homo_sapiens/GRCh38/bwa/hsa.GRCh38 \
	--no-mark-duplicates \
	--output-metadata-file=${configfile} \
	--mode=slurm \
	${fastqdir}/${id}.*.fq.gz

