#!/bin/bash
#SBATCH --partition=general
#SBATCH --mem=1024
#SBATCH --job-name publish_ge
#SBATCH --output publish.out

if [ ! $1 ]
then
    >&2 echo Project identifier is required.
    exit 1
fi

if [[ ! $1 =~ GEP[0-9]{5} ]]
then
    >&2 echo $1 does not look like a GE project identifier.
    exit 1
fi


echo 'Running copy to bioinf-ge001'

rsync -av --delete-excluded \
--exclude=temp --exclude=JobOutputs --exclude=*.out --exclude=logs \
./ ge@bioinf-ge001.cri.camres.org:/data/ge/$1/

echo 'Job done'
