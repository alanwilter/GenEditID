Running the Pipeline
====================

Fetching data
-------------

```
/home/bioinformatics/software/pipelines/kickstart/kickstart-1.5.0/bin/kickstart \
--aligner=bwamem --species=homo_sapiens --genome-version=GRCh38_hs38d1 \
--library=SLX-15103 --fastq-only --combine-samples
```

Once this has finished and the `realignment-meta.xml` pipeline driver file has been
created, the driver file must be edited to NOT mark duplicates. Add to its `<variables>`
block:

```XML
    <markDuplicates>false</markDuplicates>
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

/home/bioinformatics/software/pipelines/alignment/alignment-2.3.0/bin/run-pipeline \
--mode=slurm \
/mnt/scratchb/bioinformatics/bowers01/CRISPR/realignment-meta.xml
```

Amplicon Pipeline
-----------------

### Create driver files

Use the `create_pipeline_files.py` script to create `amplicons.txt` and `targets.txt`,
and convert the sample sheet CSV from kickstart to `samples.txt`:

```
python scripts/create_pipeline_files.py --project="GEP00004" \
--seq-dict="/mnt/scratchb/bioinformatics/reference_data/reference_genomes/homo_sapiens/GRCh38_hs38d1/fasta/hsa.GRCh38_hs38d1.dict" \
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
