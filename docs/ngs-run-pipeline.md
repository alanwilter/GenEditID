# Run the Amplicon Sequencing Pipeline

## Dependencies

The following software needs to be installed:

- [x] Unix Bourne and Bash shells
- [x] Java SE 8 or above, add it to your `.bashrc` file, and test java version `java -version`
  ```
  export JAVA_HOME=/home/bioinformatics/software/java/latest
  export PATH=${JAVA_HOME}/bin:${PATH}
  ```
- [x] Genome Analysis Toolkit (GATK) [version 3.7]
  ```
  /home/bioinformatics/software/pipelines/ampliconseq/ampliconseq-pipeline-0.4.1/gatk.jar
  ```
- FreeBayes (optional)
- [x] VarDict (Java version only, optional) [version 1.5.1], set `<vardictExecutable>` to
  ```
  /home/bioinformatics/software/vardict/VarDict-1.5.1/bin/VarDict
  ```
- [x] R including the packages Nozzle.R1, base64, scales, forcats, readr, tidyr, dplyr, ggplot2, fitdistrplus
  ```
  /home/bioinformatics/pipelinesoftware/ampliconseq/el7/R-3.4.0/lib64/R/library/
  ```
- R packages shiny, DT and highcharter for R/Shiny visualization tool (optional)
- [x] Perl
  ```
  cd ~
  wget -O - https://install.perlbrew.pl | bash
  perl5/perlbrew/bin/perlbrew install perl-5.16.0
  perl5/perlbrew/bin/perlbrew use perl-5.16.0
  curl -L https://cpanmin.us | perl - App::cpanminus
  mkdir tmp
  ```
  add these lines to your `.bashrc` file
  ```
  export PERLBREW_ROOT=${HOME}/perl5/perlbrew
  export PERLBREW_HOME=${HOME}/tmp/.perlbrew
  source ${PERLBREW_ROOT}/etc/bashrc
  ${PERLBREW_ROOT}/bin/perlbrew use perl-5.16.0
  ```
- [ ] [Ensembl Variant Effect Predictor, release 89](http://www.ensembl.org/info/docs/tools/vep/script/vep_download.html) and dependent Perl packages (File::Copy::Recursive, Archive::Zip, DBI)
  ```
  perl5/perlbrew/bin/perlbrew use perl-5.16.0
  cpanm Archive::Zip
  cpanm DBI
  cpanm File::Copy::Recursive
  cpanm Bio::Root::Version
  curl -L -O https://github.com/Ensembl/ensembl-vep/archive/release/89.zip
  unzip 89.zip
  cd ensembl-vep-release-89/
  perl INSTALL.pl  
  # ----------------------------------------------------------------------------
  # create cache files in /home/pajon01/.vep
  # specify species files 45 59
  # 45 : homo_sapiens_refseq_vep_89_GRCh38.tar.gz
  # 59 : mus_musculus_refseq_vep_89_GRCm38.tar.gz
  ```


## Reference genome

Use `/scratchb/bioinformatics/reference_data/reference_genomes/homo_sapiens/GRCh38_hs38d1`, a GATK-compatible reference genome which is the NCBI reference genome without alt sequences with the decoy and EBV. The hs38d1 part of the name reflects that this reference includes the decoy sequence that has the name hs38d1.


## Copy scripts and config files onto cluster in project folder

- Create GE project folder on cluster
```
ssh clust1-headnode
cd /path/to/scratch/space
mkdir GEPID
```

- Copy scripts and ampliconseq config files
```
scp shell/ngs/* clust1-headnode:/path/to/scratch/space/GEPID/.
```

## Fetch fastq files using kickstart

See [Kickstart documentation](https://intranet.cri.camres.org/core-facilities/bioinformatics/sequencing/kickstart)
It is intended to help with setting up a consistent working directory structure and to fetch and prepare files for downstream analysis work which generates a meta data file to use with the alignment pipeline to align FASTQ data.

```
sbatch job_kickstart.sh SLX-ID
```


## Align

- Edit `realignment-meta.xml` by adding `<markDuplicates>false</markDuplicates>` into the list of `<variables>`.

```
sbatch job_alignment.sh
```


## Run Amplicon sequencing pipeline

### Configure

Create directories for analysis, replace GEPID in config file by current one, and make the BAMs match the ids:
```
./configure_amplicon.sh GEPID
```

### Create targets, amplicons and samples files from GE database

Get dict file from the reference genome:
`/scratchb/bioinformatics/reference_data/reference_genomes/homo_sapiens/GRCh38_hs38d1/fasta/hsa.GRCh38_hs38d1.dict`

Extract amplicons and targets coordinates from the database using script `create_pipeline_files.py`, and generate samples file:
```
python python/scripts/create_pipeline_files.py --project=GEPID --seq-dict=/path/to/hsa.GRCh38.dict --samplesheet=/path/to/samplesheet.csv
```
And copy these three files into your project directory.

### Run the amplicon pipelines

```
sbatch job_amplicon_vardict.sh
sbatch job_amplicon_gatk.sh
```
