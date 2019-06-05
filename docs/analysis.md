# EditID Analysis

## Install dependencies

- [`fastq-join`](https://github.com/brwnj/fastq-join)
- [`seqkit`](https://github.com/shenwei356/seqkit)
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [MultiQC](https://multiqc.info/)
NB. These tools need to be on your path to be executable.

## Install and run the EditID tools


### Setup and configuration

```
git clone https://github.com/GenEditID/editID.git
cd editID/
python3 -m venv venv
source venv/bin/activate
pip install -r python/requirements.txt
```

- Create the database

- Start the WebApp

- Create a project in the WebApp, retrieve its GEPID and load its project layout excel file into the backend database either using the WebApp or a script e.g. `shell/load_project_GEP00001.sh`

- Create a GE project folder
```
cd /path/to/my/data/
mkdir GEPID
```

- Copy all NGS scripts onto project folder
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
./job_joinreads.sh
```
- or merged when target size is smaller than the read length (`seqkit` needs to be installed)
```
cd /path/to/my/data/GEPID/
./job_mergereads.sh
```


### Read counts

- Get fasta file from the reference genome `hsa.GRCh38_hs38d1.fa`

- Extract amplicons and targets coordinates from the database using script `create_pipeline_files.py`, and config file:
  ```
  cd /path/to/my/data/GEPID/
  source ~/editID/venv/bin/activate
  python ~/editID/python/scripts/create_ampli_count_conf.py --project=GEPID --genome=/path/to/hsa.GRCh38_hs38d1.fa
  ```

- Run [`run_ampli_count.py`](../python/scripts/run_ampli_count.py) script on all fasta files
  ```
  cd /path/to/my/data/GEPID/
  ./job_amplicount.sh
  tail -f amplicount.out
  ```

- Check results in `amplicount.csv`


### Identify variants and plot results

- Run [`run_variant_id.py`](../python/scripts/run_variant_id.py) script from the project directory:
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
./job_fastqc.sh
```

- When alignment is done, run MultiQC report:
```
cd /path/to/my/data/GEPID/
./job_multiqc.sh
```
NB. You need in your home directory a virtual environment called `venv-multiqc` with MultiQC installed to be able to run it.
