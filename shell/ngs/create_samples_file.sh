#!/bin/bash

if [ $# -lt 2 ]
then
	echo "Usage: $0 <sample_sheet_file> <samples_file>"
	exit 1
fi

samplesheet=$1
samplesfile=$2

echo "ID	SAMPLE" > $samplesfile

sed '1d;s/^"//;s/"$//;s/","/	/g' $samplesheet | awk 'BEGIN { FS = "\t"; OFS = "\t" } { print $3, $1 }' | sort >> $samplesfile

