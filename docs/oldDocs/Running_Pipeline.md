Running the Pipeline
====================

Fetching data
-------------

```
/home/bioinformatics/software/pipelines/kickstart/current/bin/kickstart \
--no-mark-duplicates --aligner=bwamem --genome=GRCh38_hs38d1 --library=SLX-15103
```

Alignment
---------

Create a script for Slurm:

```
#!/bin/bash
#SBATCH --no-requeue
#SBATCH -p general
#SBATCH -J SLX-15103
#SBATCH --nodes 1
#SBATCH --mem 4096
#SBATCH --mincpus 1
#SBATCH --open-mode truncate
#SBATCH -o alignment.log

/home/bioinformatics/software/pipelines/alignment/current/bin/run-pipeline alignment-meta.xml
```

Amplicon Pipeline
-----------------

### Create driver files

Use the `create_pipeline_files.py` script to create `amplicons.txt` and `targets.txt`,
and convert the sample sheet CSV from kickstart to `samples.txt`:

```
python scripts/create_pipeline_files.py --project="GEP00004" \
--seq-dict="/mnt/scratcha/bioinformatics/reference_data/reference_genomes/homo_sapiens/GRCh38_hs38d1/fasta/hsa.GRCh38_hs38d1.dict" \
--samplesheet="samplesheet.csv"
```

### Make the BAMs match the ids

```
mkdir idbam
cd idbam
for bam in ../bam/*.bam
do
    file=$(basename $bam)
    barcode=$(echo $file | cut -d '.' -f 2).bam
    ln -s $bam $barcode
done
for bai in ../bam/*.bai
do
    file=$(basename $bai)
    barcode=$(echo $file | cut -d '.' -f 2).bai
    ln -s $bai $barcode
done
```
