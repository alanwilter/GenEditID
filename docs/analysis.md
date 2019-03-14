# Genome Editing Analysis

## Table of contents
- [Installation instructions](#Installation-instructions)
- [Steps to setup a project and run the analysis](#Steps-to-setup-a-project-and-run-the-analysis)


## Installation instructions

### Join/merge reads

Install:
- [`fastq-join`](https://github.com/brwnj/fastq-join)
- [`seqkit`](https://github.com/shenwei356/seqkit)

### FastQC
- Installation
```
ssh clust1-headnode
cd bin
ln -s /home/bioinformatics/software/fastqc/fastqc-v0.11.8/fastqc
```

- Usage
```
scp shell/ngs/job_fastqc.sh clust1-headnode:/scratchb/bioinformatics/pajon01/genome-editing/GEP00010/.
ssh clust1-headnode
cd /scratchb/bioinformatics/pajon01/genome-editing/GEP00010/
sbatch job_fastqc.sh
```

### MultiQC
https://multiqc.info/

- Installation
```
ssh clust1-headnode
/home/bioinformatics/software/python/python-3.4.9/bin/python3 -m venv venv-multiqc
source venv-multiqc/bin/activate
pip install multiqc
```

- Usage
```
scp shell/ngs/job_multiqc.sh clust1-headnode:/scratchb/bioinformatics/pajon01/genome-editing/GEP00010/.
ssh clust1-headnode
cd /scratchb/bioinformatics/pajon01/genome-editing/GEP00010/
sbatch job_multiqc.sh
```

### CRISPResso
https://github.com/lucapinello/CRISPResso

- Installing
```
```

- Usage
```
```

### AmpliCan

http://bioconductor.org/packages/release/bioc/html/amplican.html

- Installation
```
```

- Usage
```
```

### Find amplicon coordinates

Install:
- `faToTwoBit`
- `gfServer` from the BLAT tool suite

Investigate: https://github.com/mdshw5/pyfaidx

### Install the Amplicon Sequencing Pipeline

- Pipeline repository: It is located under the Bioinformatics Core subversion repository in
`svn+ssh://svn-bioinformatics.cruk.cam.ac.uk/data/mib-cri/SVNREP/pipelines/ampliconseq/trunk`

- Dependencies have been install for the CRUK-CI cluster, and detailed information can be found in the `docs/INSTALL.md` pipeline repository. These are located into `/home/bioinformatics/pipelinesoftware/ampliconseq/el7/`
  - The following software needs to be installed:
    - Unix Bourne and Bash shells
    - Java SE 8 or above, add it to your `.bashrc` file, and test java version `java -version`
      ```
      export JAVA_HOME=/home/bioinformatics/software/java/latest
      export PATH=${JAVA_HOME}/bin:${PATH}
      ```
    - Genome Analysis Toolkit (GATK) [version 3.7]
      ```
      /home/bioinformatics/software/pipelines/ampliconseq/ampliconseq-pipeline-0.4.1/gatk.jar
      ```
    - FreeBayes (optional)
    - VarDict (Java version only, optional) [version 1.5.1], set `<vardictExecutable>` to
      ```
      /home/bioinformatics/pipelinesoftware/ampliconseq/el7/VarDict-1.5.1/bin/VarDict
      ```
    - R including the packages Nozzle.R1, base64, scales, forcats, readr, tidyr, dplyr, ggplot2, fitdistrplus
      ```
      /home/bioinformatics/pipelinesoftware/ampliconseq/el7/R-3.4.0/lib64/R/library/
      ```
    - R packages shiny, DT and highcharter for R/Shiny visualization tool (optional)
    - Perl
      - if you wish to install your own:
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
      - or use the one already installed, and add these lines to your `.bashrc` file
        ```
        export SOFTWARE_ROOT=/home/bioinformatics/pipelinesoftware/ampliconseq/el7
        export PERLBREW_ROOT=${SOFTWARE_ROOT}/perlbrew
        source ${PERLBREW_ROOT}/etc/bashrc
        ${PERLBREW_ROOT}/bin/perlbrew use perl-5.16.0  
        ```
    - [Ensembl Variant Effect Predictor, release 89](http://www.ensembl.org/info/docs/tools/vep/script/vep_download.html) and dependent Perl packages (File::Copy::Recursive, Archive::Zip, DBI)
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
      # specify species files
      # 47 : homo_sapiens_vep_89_GRCh38.tar.gz
      # 60 : mus_musculus_vep_89_GRCm38.tar.gz
      # FASTA files for the following species are available
      # 29 : homo_sapiens
      # 38 : mus_musculus
      ```
      and add `vep` to your path `ln -s /home/pajon01/ensembl-vep-release-89/vep ~/bin/.`
  - In summary, if using the CRUKCI infrastructure, you only need to add these lines to your `.bashrc`:
    ```
    ### AmpliconSeq pipeline dependencies

    export AMPLICONSEQ_SOFT_ROOT=/home/bioinformatics/pipelinesoftware/ampliconseq/el7

    # java
    export JAVA_HOME=/home/bioinformatics/software/java/latest
    export PATH=${JAVA_HOME}/bin:${PATH}

    # perl
    export PATH=${AMPLICONSEQ_SOFT_ROOT}/perl-5.16.0/bin:${PATH}

    # ensembl vep
    export VEP_CACHE=/mnt/scratchb/bioinformatics/reference_data/ensembl/vep
    export PATH=${AMPLICONSEQ_SOFT_ROOT}/ensembl-vep-release-91.1/htslib:${PATH}
    export PATH=${AMPLICONSEQ_SOFT_ROOT}/ensembl-vep-release-91.1/:${PATH}
    ```

- Reference genome: use `/scratcha/bioinformatics/reference_data/reference_genomes/homo_sapiens/GRCh38_hs38d1`, a GATK-compatible reference genome which is the NCBI reference genome without alt sequences with the decoy and EBV. The hs38d1 part of the name reflects that this reference includes the decoy sequence that has the name hs38d1.


## Steps to setup a project and run the analysis


### Step1: Setup and configuration

- Create a project in webapp, retrieve its GEPID and load its project layout excel file into the backend database either using the webapp or a script e.g. `shell/load_project_GEP00013.sh`

- Submit project for sequencing using the webapp which will automatically add the samples into the acceptance queue in Clarity. Retrieve project name and SLXID.

- Create Redmine project and copy id into Clarity to be able to receive notification email when sequencing data is available

- Create a GE project folder on cluster
```
ssh clust1-headnode
cd /path/to/scratch/space
mkdir GEPID
```

- Copy all NGS scripts onto project folder on cluster
```
scp shell/ngs/* clust1-headnode:/path/to/scratch/space/GEPID/.
```


### Step2: Fetch fastq files

See [Kickstart documentation](http://internal-bioinformatics.cruk.cam.ac.uk/docs/kickstart/usage.html)
It is intended to help with setting up a consistent working directory structure and to fetch and prepare files for downstream analysis work which generates a meta data file to use with the alignment pipeline to align FASTQ data.

#### 2.1 Download files
```
sbatch job_kickstart.sh SLX-ID
```

#### 2.2 Check sequencing information
- SLXID and
- Read length PE150 or PE300

#### 2.3 Combine reads

- Reads should be joined when target size is bigger than read length (`fastq-join` needs to be installed)
```
sbatch job_joinreads.sh
```
- or merged when target size is smaller than the read length (`seqkit` needs to be installed)
```
sbatch job_mergereads.sh
```


### Step3: Read counts



### Step4: Align

```
sbatch job_alignment.sh
```


### Step5: MultiQC report

- Run FastQC on all joined fastq sample files:
```
sbatch job_fastqc.sh
```

- When alignment is done, run MultiQC report:
```
sbatch job_multiqc.sh
```
NB. You need in your home directory on the cluster a virtual environment with MultiQC installed to be able to run it.


### Step6: Run Amplicon sequencing pipeline when alignment done

#### 6.1 Configure

Create directories for analysis, replace GEPID in config file by current one, and make the BAMs match the ids:
```
./configure_amplicon.sh GEPID
```

#### 6.2 Create targets, amplicons and samples files from GE database

Get dict file from the reference genome:
`/scratcha/bioinformatics/reference_data/reference_genomes/homo_sapiens/GRCh38_hs38d1/fasta/hsa.GRCh38_hs38d1.dict`

Extract amplicons and targets coordinates from the database using script `create_pipeline_files.py`, and generate samples file:
```
python python/scripts/create_pipeline_files.py --project=GEPID --seq-dict=/path/to/hsa.GRCh38.dict --filelist=/path/to/filelist.csv
```
And copy these three files (`targets.txt`, `amplicons.txt` and `samples.txt`) into your project directory.

#### 6.3 Run the amplicon pipelines

```
sbatch job_amplicon_vardict.sh
sbatch job_amplicon_gatk.sh
```

Check the output folders for results in `amplicon_gatk/output`:
- First, look at the coverage report `*.amplicon_coverage.xlsx` to see if reads have been aligned to the targets
- Then, check variants in `*.variants.xlsx`
- If data looks alright, it can then be loaded into the database


### Step7: Load results into database and calculate score

#### 7.1 Load pipeline variant result files

```
python python/scripts/load_variant_results.py --file=amplicon_gatk/output/GEP00012.gatk.variants.xlsx --project_geid=GEP00012 --caller=HaplotypeCaller
python python/scripts/load_variant_results.py --file=amplicon_vardict/output/GEP00012.vardict.variants.xlsx --project_geid=GEP00012 --caller=VarDict
```

#### 7.2 Calculate score

```
python python/scripts/load_mutation_summary.py --project_geid=GEP00012
```


### Step8: Copy and Backup

- Copy data to the GE server

```
sbatch job_publish.sh <project_geid>
```

This copies everything from the processing directory except working (temporary)
directories and log files or directories.

- Give access to BAM files and create links in result table
- Backup analysis files
- Archive data
