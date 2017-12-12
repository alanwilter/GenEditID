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

Note reference genomes in the standard path `/scratchb/bioinformatics/reference_data/reference_genomes/homo_sapiens/GRCh38` are not GATK compatible, therefore create your own to be compatible with GATK in `/mnt/scratcha/bioinformatics/chilam01/reference_data/homo_sapiens/GRCh38`.

```
sh create_alignment_config.sh SLX-13775 /mnt/scratcha/bioinformatics/chilam01/projects/CRISPR/NGS/SLX-13775/run_and_document/fastq  general alignment.xml
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
For instance for GRCh38 header is from `/scratchb/bioinformatics/reference_data/reference_genomes/homo_sapiens/GRCh38/fasta/hsa.GRCh38.dict`

Extract amplicon and target coordinates from the database using script `???`

### Create sample name and id mapping file

```
sh create_samples_file.sh samplesheet.csv samples.txt
```

### Haplotype caller config

- Copy haplotype config file and modify
```
cp /home/bioinformatics/software/pipelines/ampliconseq/ampliconseq-pipeline-0.4/config/config.gatk.xml .
```

- Edit these lines

```
<referenceDataDir>/mnt/scratcha/bioinformatics/chilam01/projects/CRISPR/NGS/SLX-13775/test_new_pip_version_and_GRCh38/refDir/GRCh38</referenceDataDir>

<ampliconIntervals>${referenceDataDir}/amplicons.txt</ampliconIntervals>
<targetIntervals>${referenceDataDir}/targets.txt</targetIntervals>

<referenceSequence>/mnt/scratcha/bioinformatics/chilam01/reference_data/homo_sapiens/GRCh38/hsa.GRCh38.fa</referenceSequence>
```

### Run the amplicon pipeline

```
# perlbrew was used for bioperl instalation
# therefore run this command before you start the amplicon pipeline
perlbrew use perl-5.16.0

sbatch job_amplicon_pipeline.sh
```
