# Run the Amplicon Sequencing Pipeline


## Fetch fastq files

```
sbatch job_kickstart.sh
```


## Alignment steps

```
mkdir bam
cd bam
```

### Create alignment configuration file

Note reference genomes in the standard path `/scratchb/bioinformatics/reference_data/reference_genomes/homo_sapiens/GRCh38` are not GATK compatible, therefore use Matt's GATK compatible ones in `/scratchb/bioinformatics/eldrid01/reference_data/reference_genomes/homo_sapiens/GRCh38/`.

```
./create_alignment_config.sh SLX-13775 /path/to/project/fastq  general alignment.xml
```

Modify `alignment.xml`:
- remove line with samtools
- add software path `<softwareDir>/home/bioinformatics/pipelinesoftware/alignment/el7</softwareDir>`

### Run alignment pipeline

```
sbatch job_alignment.sh
```


## Amplicon sequencing steps

### Create target and amplicon files

Both amplicon and target files should be in picard style with header section
Get header from .dict file from reference genomes
For instance for GRCh38 header is from `/scratchb/bioinformatics/eldrid01/reference_data/reference_genomes/homo_sapiens/GRCh38/fasta/hsa.GRCh38.dict`

Extract amplicon and target coordinates from the database using script `create_pipeline_files.py`:
```
python python/scripts/create_pipeline_files.py --project=GEP00006 --seq-dict=data/seqdict/hsa.GRCh38.dict
```

And copy these files into your project directory under new folder `amplicon`.

### Create sample name and id mapping file

```
./create_samples_file.sh samplesheet.csv samples.txt
```

### Haplotype caller config

- Copy haplotype config file and modify
```
cp /home/bioinformatics/software/pipelines/ampliconseq/ampliconseq-pipeline-0.4.1/config/config.gatk.xml /path/to/project/amplicon/.
```

- Edit these lines

```
<mode>slurm</mode>
<referenceSequence>/scratchb/bioinformatics/eldrid01/reference_data/reference_genomes/homo_sapiens/GRCh38/fasta/hsa.GRCh38.fa</referenceSequence>
<variantEffectPredictorAssembly>GRCh38</variantEffectPredictorAssembly>
```

### Run the amplicon pipeline

- Install `perl` dependency using `perlbrew` prior to run the pipeline

```
cd ~
curl -L https://install.perlbrew.pl | bash
ln -s ~/perl5/perlbrew/bin/perlbrew bin/.
perlbrew install perl-5.16.0
```

- Run the pipeline
```
sbatch job_amplicon_pipeline.sh
```
