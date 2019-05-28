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

### Find amplicon coordinates

Install:
- `faToTwoBit`
- `gfServer` from the BLAT tool suite

Investigate: https://github.com/mdshw5/pyfaidx

- Reference genome: use `/scratcha/bioinformatics/reference_data/reference_genomes/homo_sapiens/GRCh38_hs38d1`, a GATK-compatible reference genome which is the NCBI reference genome without alt sequences with the decoy and EBV. The hs38d1 part of the name reflects that this reference includes the decoy sequence that has the name hs38d1.


### EditID tool

```
git clone https://github.com/GenEditID/editID.git
cd editID/
python3 -m venv venv
source venv/bin/activate
pip install -r python/requirements.txt
```

## Steps to setup a project and run the analysis


### Setup and configuration

- Create a project in webapp, retrieve its GEPID and load its project layout excel file into the backend database either using the webapp or a script e.g. `shell/load_project_GEP00001.sh`

- Create a GE project folder on cluster
```
ssh clust1-headnode
cd /path/to/my/data/
mkdir GEPID
```

- Copy all NGS scripts onto project folder on cluster
```
cp ~/editID/shell/ngs/* /path/to/my/data/GEPID/.
```


### Fetch fastq files

#### Check sequencing information
- Read length PE150 or PE300

#### Combine reads

- Reads should be joined when target size is bigger than read length (`fastq-join` needs to be installed)
```
cd /path/to/my/data/GEPID/
sbatch job_joinreads.sh
```
- or merged when target size is smaller than the read length (`seqkit` needs to be installed)
```
cd /path/to/my/data/GEPID/
sbatch job_mergereads.sh
```


### Read counts

- Get fasta file from the reference genome:
  - `/scratcha/bioinformatics/reference_data/reference_genomes/homo_sapiens/GRCh38_hs38d1/fasta/hsa.GRCh38_hs38d1.fa`

- Extract amplicons and targets coordinates from the database using script `create_pipeline_files.py`, and config file:
  ```
  cd /path/to/my/data/GEPID/
  source ~/editID/venv/bin/activate
  python ~/editID/python/scripts/create_ampli_count_conf.py --project=GEPID --genome=/path/to/hsa.GRCh38_hs38d1.fa
  ```


- Run amplicount script on all fasta files
  ```
  cd /path/to/my/data/GEPID/
  sbatch job_amplicount.sh
  tail -f amplicount.out
  ```

- Check results in `amplicount.csv`

- Classify variants and plot results from the project directory:
  ```
  cd /path/to/my/data/GEPID/
  source ~/editID/venv/bin/activate
  python ~/editID/python/scripts/run_variant_id.py
  ```

- Check results in `editid_variantid/variantid.csv` and `editid_variantid/impacts.csv` and plots
  - `editid_variantid/coverage_[amplicon_id].html`
  - `editid_variantid/koscores_[amplicon_id].html`

- Retrieve sample location on plates from the database and add them onto the necessary files
```
cd /path/to/my/data/GEPID/
source ~/editID/venv/bin/activate
python ~/editID/python/scripts/get_sample_loc.py GEPID
python ~/editID/python/scripts/add_sample_loc.py
```

- Plot heatmap on plates
```
cd /path/to/my/data/GEPID/
source ~/editID/venv/bin/activate
python ~/editID/python/scripts/plot_scores.py
```

- Visualise plots
  - `editid_variantid/heatmap_[amplicon_id].html`
  - `editid_variantid/heatmap_protein_expression.html` (if data available)
  - `editid_variantid/heatmap_combined_data.html` (if data available)


### MultiQC report

- Run FastQC on all joined fastq sample files:
```
cd /path/to/my/data/GEPID/
sbatch job_fastqc.sh
```

- When alignment is done, run MultiQC report:
```
cd /path/to/my/data/GEPID/
sbatch job_multiqc.sh
```
NB. You need in your home directory on the cluster a virtual environment with MultiQC installed to be able to run it.


### Load results into database

#### Load results
