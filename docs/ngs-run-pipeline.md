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

### Create target and amplicon files from GE database

Get header from .dict file from the reference genome:
`/scratchb/bioinformatics/reference_data/reference_genomes/homo_sapiens/GRCh38_hs38d1/fasta/hsa.GRCh38_hs38d1.dict`

Extract amplicon and target coordinates from the database using script `create_pipeline_files.py`:
```
python python/scripts/create_pipeline_files.py --project=GEP00006 --seq-dict=/path/to/hsa.GRCh38.dict
```

And copy these files into your project directory `amplicon/` folder.

### Create sample name and id mapping file

```
./create_samples_file.sh samplesheet.csv amplicon/samples.txt
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


```
mkdir amplicon/vardict
mkdir amplicon/gatk
sbatch job_amplicon_vardict.sh
sbatch job_amplicon_gatk.sh
```
