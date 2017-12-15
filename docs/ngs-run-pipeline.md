# Run the Amplicon Sequencing Pipeline

## Reference genome

Use `/scratchb/bioinformatics/reference_data/reference_genomes/homo_sapiens/GRCh38_hs38d1`, a GATK-compatible reference genome which is the NCBI reference genome without alt sequences with the decoy and EBV. The hs38d1 part of the name reflects that this reference includes the decoy sequence that has the name hs38d1.


## Copy scripts onto cluster in project folder

- Create GE project folder on cluster
```
ssh clust1-headnode
cd /path/to/scratch/space
mkdir GEP00006
```

- Copy scripts and ampliconseq config files
```
scp shell/ngs/* clust1-headnode:/path/to/scratch/space/GEP00006/.
```

## Fetch fastq files using kickstart

See [Kickstart documentation](https://intranet.cri.camres.org/core-facilities/bioinformatics/sequencing/kickstart)
It is intended to help with setting up a consistent working directory structure and to fetch and prepare files for downstream analysis work which generates a meta data file to use with the alignment pipeline to align FASTQ data.

Edit `job_kickstart.sh` to update `--library=SLX-ID` to the one associated with the GE project
```
sbatch job_kickstart.sh
```


## Alignment

```
sbatch job_alignment.sh
```


## Amplicon sequencing steps

### Create targets, amplicons and samples files from GE database

Get dict file from the reference genome:
`/scratchb/bioinformatics/reference_data/reference_genomes/homo_sapiens/GRCh38_hs38d1/fasta/hsa.GRCh38_hs38d1.dict`

Extract amplicons and targets coordinates from the database using script `create_pipeline_files.py`, and generate samples file:
```
python python/scripts/create_pipeline_files.py --project=GEP00006 --seq-dict=/path/to/hsa.GRCh38.dict --samplesheet=/path/to/samplesheet.csv
```

And copy these files into your project directory `amplicon/` folder.

### Make the BAMs match the ids

```
./rename_bam_files.sh
```

### Vardict config

Edit this line in `config.vardict.xml`

```
<runId>GEPID</runId>
```

### Gatk config

Edit this line in `config.gatk.xml`

```
<runId>GEPID</runId>
```


### Run the amplicon pipelines

NB. You need to have Java 8 installed and you may have to add it to your path.
```
export JAVA_HOME=/home/bioinformatics/software/java/latest
export PATH=${JAVA_HOME}/bin:${PATH}
```

```
mkdir amplicon/vardict
mkdir amplicon/gatk
sbatch job_amplicon_vardict.sh
sbatch job_amplicon_gatk.sh
```
