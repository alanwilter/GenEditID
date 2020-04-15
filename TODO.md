# Things done...

- FastQ files added into GDrive to be [downloaded from here](https://drive.google.com/drive/folders/1MN_vzy3hjGOAnycwwtI53nrAuOaB5RJf?usp=sharing) (step 3)
- Access to [Google Drive paper from here](https://drive.google.com/drive/folders/1MQAmhxjuewH2gDoUkzXzz1wmgMK6CV7E?usp=sharing)
- Slack integration, find us on [geneditid](geneditid.slack.com)

# Things in progress...

- Create proper test data
  - [x] select 3 FastQ files
  - [x] run all analysis steps on these 3 files only
  - [ ] improve setup for testing

- Add new plots to the WebApp
  - [ ] Update WebApp interface to reflect changes
  - [ ] Add output of the run_ampli_count tool into a new database table for plotting (add class AmpliCountResult in model.py)
  - [ ] Update the plotting scripts to call the code from the WebApp to avoid code duplication


# Things to do...

- Simplify submission spreadsheet and database

- Add validation scripts for the parameters chosen for the alignment `pairwise2.align.globalms(ref_sequence, variant['sequence'], 5, -4, -3, -0.1)`
- Add validation scripts for variant classification `Variant(contig=contig, start=start, ref=ref, alt=alt, ensembl=genome)` `var.effects().top_priority_effect()`
- Modify the analysis steps to be ran from the database instead of flat files
- Drive the analysis from the WebApp

# Slack integration

- Find us on [slack](geneditid.slack.com)

# Test data: select 3 FastQ files

The test data is from the project GEP00009, it includes fastq run files for FLD0018, FLD0046, FLD0060.
The target is FTO on chr16 chr16_53704130, length 182.
In folder [test_data/GEP00009/fastq/](https://drive.google.com/drive/folders/1TCYqPAkrP6ju-lop_l2Ymb4lSsb1oHtq?usp=sharing) in GDrive, there are the 6 files needed for testing :
- SLX-15026.FLD0018.000000000-BJWVR.s_1.r_1.fq.gz
- SLX-15026.FLD0018.000000000-BJWVR.s_1.r_2.fq.gz
- SLX-15026.FLD0046.000000000-BJWVR.s_1.r_1.fq.gz
- SLX-15026.FLD0046.000000000-BJWVR.s_1.r_2.fq.gz
- SLX-15026.FLD0060.000000000-BJWVR.s_1.r_1.fq.gz
- SLX-15026.FLD0060.000000000-BJWVR.s_1.r_2.fq.gz
These files are paired-end reads of length 300.

# Run all analysis steps

From [step 3](https://geneditid.github.io/manual.html) onward.

## Step 3: Fetch fastq files

- Get all fastq files related to your project into the GEPID project folder into a `fastq` folder
Example of folder structure:
  ```
  cd
  mkdir -p ~/data/geneditid_testdata/GEP00009/fastq
  cd ~/data/geneditid_testdata/GEP00009/fastq
  ```
  Add these six files into this folder:
  - SLX-15026.FLD0018.000000000-BJWVR.s_1.r_1.fq.gz
  - SLX-15026.FLD0018.000000000-BJWVR.s_1.r_2.fq.gz
  - SLX-15026.FLD0046.000000000-BJWVR.s_1.r_1.fq.gz
  - SLX-15026.FLD0046.000000000-BJWVR.s_1.r_2.fq.gz
  - SLX-15026.FLD0060.000000000-BJWVR.s_1.r_1.fq.gz
  - SLX-15026.FLD0060.000000000-BJWVR.s_1.r_2.fq.gz

- Combine paired-end reads by merging or joining reads to generate `.fqjoin.gz` files for ampli_count analysis:
Here we have a target of 180 bases, and paired-end read of 300 in length, so we can merged the reads because the target size is smaller than the read length.
  ```
  cd ~/data/geneditid_testdata/GEP00001/
  ~/GenEditID/shell/ngs/job_mergereads.sh
  ```
  NB. You need `seqkit` on your path to be able to run this step. [Download it](https://bioinf.shenwei.me/seqkit/download/)!

All the following steps should be ran in the **GEPID project folder**, after activating the virtual python environment where all needed libraries have been installed:
```
cd ~/data/geneditid_testdata/GEP00001/
source ~/GenEditID/venv/bin/activate
```

## Step 4: Run `ampli_count`

- Download Reference Genome (release 95): [Homo_sapiens.GRCh38.dna.toplevel.fa.gz](ftp://ftp.ensembl.org/pub/release-95/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz) check your disk space, before unzipping the reference genome
  ```
  cd ~/data/geneditid_testdata/GEP00001/
  gunzip Homo_sapiens.GRCh38.dna.toplevel.fa.gz
  ```

- Extract amplicons and targets coordinates from the database to produce the config file to run `ampli_count`:
  NB. We are testing with data files from project GEP00009, so we need to create 8 'fake' projects first in the webapp  and then upload the submission sheet of project GEP00009 (re-do steps 1 and 2) to be able to run the script and get sensible data:
  ```
  geneditid_create_amplicount_config --project=GEP00001 --genome=Homo_sapiens.GRCh38.dna.toplevel.fa
  --- Amplicon #1
  Target name	FTO
  Forward primer	TCCAGGGCGAGGGATCTAC
  Reverse primer	GAAGGGCTTGGTTTGATGGC
  target coord	chr16	53704149	53704290	+	chr16_53704130
  amplicon coord	chr16	53704130	53704310	+	chr16_53704130
  amplicon desc	chr16:53704130:53704310
  amplicon seq	TCCAGGGCGAGGGATCTACGCAGCTTGCGGTGGCGAAGGCGGCTTTAGTGGCAGCATGAAGCGCACCCCGACTGCCGAGGAACGAGAGCGCGAAGCTAAGGTATGTCGGGCTCCCGGGGCCTGGAGATCTTCGTGCGCTGTGAGCAAGGATCAGGGAACCGGAAGGGCTTGGTTTGATGGC
  amplicon seq len	181
  >>> Reverse primer sequence different than the one submitted! [submitted: GCCATCAAACCAAGCCCTTC, found: GAAGGGCTTGGTTTGATGGC]
  ---
  --- Amplicon #2
  Target name	IRX5
  Forward primer	CCATGCCCGTGTGTG
  Reverse primer	GCTACAACTCGCACCTCCA
  target coord	chr16	54931196	54931376	+	chr16_54931181
  amplicon coord	chr16	54931181	54931395	+	chr16_54931181
  amplicon desc	chr16:54931181:54931395
  amplicon seq	CCATGCCCGTGTGTGGCCATGTCCTATCCGCAGGGCTACTTGTACCAGCCGTCCGCCTCGCTGGCGCTCTACTCGTGCCCGGCGTACAGCACCAGCGTCATTTCGGGGCCCCGCACGGATGAGCTCGGCCGCTCTTCTTCGGGCTCCGCGTTCTCGCCCTACGCTGGCTCGACTGCCTTCACGGCGCCCTCGCCGGGCTACAACTCGCACCTCCA
  amplicon seq len	215
  >>> Reverse primer sequence different than the one submitted! [submitted: TGGAGGTGCGAGTTGTAGC, found: GCTACAACTCGCACCTCCA]
  ---
  ```

  ```
  cat amplicount_config.csv
  id,fprimer,rprimer,amplicon,coord
  chr16_53704130,TCCAGGGCGAGGGATCTAC,GAAGGGCTTGGTTTGATGGC,TCCAGGGCGAGGGATCTACGCAGCTTGCGGTGGCGAAGGCGGCTTTAGTGGCAGCATGAAGCGCACCCCGACTGCCGAGGAACGAGAGCGCGAAGCTAAGGTATGTCGGGCTCCCGGGGCCTGGAGATCTTCGTGCGCTGTGAGCAAGGATCAGGGAACCGGAAGGGCTTGGTTTGATGGC,chr16:53704130-53704310
  chr16_54931181,CCATGCCCGTGTGTG,GCTACAACTCGCACCTCCA,CCATGCCCGTGTGTGGCCATGTCCTATCCGCAGGGCTACTTGTACCAGCCGTCCGCCTCGCTGGCGCTCTACTCGTGCCCGGCGTACAGCACCAGCGTCATTTCGGGGCCCCGCACGGATGAGCTCGGCCGCTCTTCTTCGGGCTCCGCGTTCTCGCCCTACGCTGGCTCGACTGCCTTCACGGCGCCCTCGCCGGGCTACAACTCGCACCTCCA,chr16:54931181-54931395
  ```
- Run [geneditid_run_amplicount](https://github.com/GenEditID/GenEditID/blob/master/python/scripts/run_ampli_count.py) script on all fasta files from your project directory
  ```
  geneditid_run_amplicount --fastqdir=fastq/
  ```
  Three files analysed in one hour:
  ```
  2020-04-02 23:15:10,526 dnascissors              INFO    : Counting reads in file SLX-15026.FLD0018.000000000-BJWVR.s_1.fqjoin.gz
  2020-04-02 23:27:24,464 dnascissors              INFO    : FLD0018,255430,39,96494
  2020-04-02 23:27:24,758 dnascissors              INFO    : Counting reads in file SLX-15026.FLD0046.000000000-BJWVR.s_1.fqjoin.gz
  2020-04-02 23:48:40,345 dnascissors              INFO    : FLD0046,366226,49,118951
  2020-04-02 23:48:40,571 dnascissors              INFO    : Counting reads in file SLX-15026.FLD0060.000000000-BJWVR.s_1.fqjoin.gz
  2020-04-03 00:09:30,860 dnascissors              INFO    : FLD0060,313702,52,104464
  ```

- Check results in `amplicount.csv`
  ```
  cat amplicount.csv
  sample_id,amplicon_id,total_reads,amplicon_reads,amplicon_filtered_reads,amplicon_low_quality_reads,amplicon_primer_dimer_reads,amplicon_low_abundance_reads,variant_reads,variant_frequency,sequence
  FLD0018,chr16_54931181,255430,170277,96494,73610,173,18932,47091,48.80,CCATGCCCGTGTGTGGCCATGTCCTATCCGCAGGGCTACTTGTACCAGCCGTCCGCCTCGCTTGGCGCTCTACTCGTGCCCGGCGTACAGCACCAGCGTCATTTCGGGGCCCCGCACGGATGAGCTCGGCCGCTCTTCTTCGGGCTCCGCGTTCTCGCCCTACGCTGGCTCGACTGCCTTCACGGCGCCCTCGCCGGGCTACAACTCGCACCTCCA
  FLD0018,chr16_54931181,255430,170277,96494,73610,173,18932,17682,18.32,CCATGCCCGTGTGTGCCAATTATGGGACTTACCCACCCAGATTTAGACATAGGTCAGTGGAACTGACCCTAAGAAGAGGCAGCAATATAGGTAAGAATGAAAGCTAAGGCACATCTAACAGCCATCCATGGGTGGGGAGGGCTACAACTCGCACCTCCA
  FLD0018,chr16_54931181,255430,170277,96494,73610,173,18932,91,0.09,CCATGCCCGTGTGTGCCAATTATGGGACTTACCCACCCAGATTTAGACATAGGTCAGTGGAACTGACCCTAAGAAGAGGCAGCAATATAGATAAGAATGAAAGCTAAGGCACATCTAACAGCCATCCATGGGTGGGGAGGGCTACAACTCGCACCTCCA
  FLD0018,chr16_54931181,255430,170277,96494,73610,173,18932,68,0.07,CCATGCCCGTGTGTGGCCATGTCCCATCCGCAGGGCTACTTGTACCAGCCGTCCGCCTCGCTTGGCGCTCTACTCGTGCCCGGCGTACAGCACCAGCGTCATTTCGGGGCCCCGCACGGATGAGCTCGGCCGCTCTTCTTCGGGCTCCGCGTTCTCGCCCTACGCTGGCTCGACTGCCTTCACGGCGCCCTCGCCGGGCTACAACTCGCACCTCCA
  ...
  ```

## Step 5: Identify variants and plot results

- Run [geneditid_run_variantid](https://github.com/GenEditID/GenEditID/blob/master/python/scripts/run_variant_id.py) script from your project directory:
```
geneditid_run_variantid
INFO:pyensembl.sequence_data:Loaded sequence dictionary from /Users/pajanne/Library/Caches/pyensembl/GRCh38/ensembl95/Homo_sapiens.GRCh38.cdna.all.fa.gz.pickle
INFO:pyensembl.sequence_data:Loaded sequence dictionary from /Users/pajanne/Library/Caches/pyensembl/GRCh38/ensembl95/Homo_sapiens.GRCh38.ncrna.fa.gz.pickle
INFO:pyensembl.sequence_data:Loaded sequence dictionary from /Users/pajanne/Library/Caches/pyensembl/GRCh38/ensembl95/Homo_sapiens.GRCh38.pep.all.fa.gz.pickle
Creating Amplicon Read Coverage plot for chr16_54931181
```
- Check results in `editid_variantid/variantid.csv` and `editid_variantid/impacts.csv` and plots
  - `editid_variantid/coverage_chr16_54931181.html`
  - `editid_variantid/koscores_chr16_54931181.html`

  ```
  cat editid_variantid/variantid.csv
  sample_id,amplicon_id,total_reads,amplicon_reads,amplicon_filtered_reads,amplicon_low_quality_reads,amplicon_primer_dimer_reads,amplicon_low_abundance_reads,variant_reads,variant_frequency,sequence,variant_id,variant_type,variant_consequence,variant_score
  FLD0018,chr16_54931181,255430,170277,96494,73610,173,18932,47091,48.8,CCATGCCCGTGTGTGGCCATGTCCTATCCGCAGGGCTACTTGTACCAGCCGTCCGCCTCGCTTGGCGCTCTACTCGTGCCCGGCGTACAGCACCAGCGTCATTTCGGGGCCCCGCACGGATGAGCTCGGCCGCTCTTCTTCGGGCTCCGCGTTCTCGCCCTACGCTGGCTCGACTGCCTTCACGGCGCCCTCGCCGGGCTACAACTCGCACCTCCA,var1,Insertion,FrameShift,48.8
  FLD0018,chr16_54931181,255430,170277,96494,73610,173,18932,17682,18.32,CCATGCCCGTGTGTGCCAATTATGGGACTTACCCACCCAGATTTAGACATAGGTCAGTGGAACTGACCCTAAGAAGAGGCAGCAATATAGGTAAGAATGAAAGCTAAGGCACATCTAACAGCCATCCATGGGTGGGGAGGGCTACAACTCGCACCTCCA,var2,Insertion_Deletion_Mismatch,ComplexFrameShift,16.488
  FLD0046,chr16_54931181,366226,240180,118951,121206,23,15945,35260,29.64,CCATGCCCGTGTGTGGCCATGTCCTATCCGCAGGGCTACTTGTACCAGCCGTCCGCCTCGCTTGGCGCTCTACTCGTGCCCGGCGTACAGCACCAGCGTCATTTCGGGGCCCCGCACGGATGAGCTCGGCCGCTCTTCTTCGGGCTCCGCGTTCTCGCCCTACGCTGGCTCGACTGCCTTCACGGCGCCCTCGCCGGGCTACAACTCGCACCTCCA,var1,Insertion,FrameShift,29.64
  FLD0046,chr16_54931181,366226,240180,118951,121206,23,15945,56378,47.4,CCATGCCCGTGTGTGGCCATGTCCTATCCGCAGGGCTACTTGTACCAGCCGTCCGCCTCGCTCTACTCGTGCCCGGCGTACAGCACCAGCGTCATTTCGGGGCCCCGCACGGATGAGCTCGGCCGCTCTTCTTCGGGCTCCGCGTTCTCGCCCTACGCTGGCTCGACTGCCTTCACGGCGCCCTCGCCGGGCTACAACTCGCACCTCCA,var3,Deletion,Deletion,0.0
  FLD0060,chr16_54931181,313702,208473,104464,104001,8,9652,80867,77.41,CCATGCCCGTGTGTGGCCATGTCCTATCCGCAGGGCTACTTGTACCAGCCGTCCGCCTCGCTTGGCGCTCTACTCGTGCCCGGCGTACAGCACCAGCGTCATTTCGGGGCCCCGCACGGATGAGCTCGGCCGCTCTTCTTCGGGCTCCGCGTTCTCGCCCTACGCTGGCTCGACTGCCTTCACGGCGCCCTCGCCGGGCTACAACTCGCACCTCCA,var1,Insertion,FrameShift,77.41
  ```

  ```
  cat editid_variantid/impacts.csv
  sample_id,amplicon_id,impact,impact_frequency
  FLD0018,chr16_54931181,HighImpact,48.8
  FLD0018,chr16_54931181,MediumImpact,18.32
  FLD0046,chr16_54931181,HighImpact,29.64
  FLD0046,chr16_54931181,LowImpact,47.4
  FLD0060,chr16_54931181,HighImpact,77.41
  ```

## Step 6: Plot koscore heatmap on plates

- Retrieve sample location on plates from the database and add them onto the necessary files
  ```
  geneditid_add_sample_location GEPID
  ```

- Plot heatmap on plates
  ```
  geneditid_plot_scores
  ```

- Visualise plots in your browser
  - `editid_variantid/heatmap_chr16_54931181.html`


# Issue with large uncompressed genome file needed to identify coordinates of amplicon

- Get [Homo_sapiens.GRCh38.dna.toplevel.fa.gz](ftp://ftp.ensembl.org/pub/release-95/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz)
- [Read code from pyfaidx](https://github.com/mdshw5/pyfaidx/blob/master/tests/test_Fasta_bgzip.py): it seems possible to give a gzip file to Fasta() method
- Try directly but without success
  ```
  geneditid_create_amplicount_config --project=GEP00009 --genome=/Users/pajanne/workspace/data/Homo_sapiens.GRCh38.dna.toplevel.fa.gz
  --- Amplicon #1
  Target name	FTO
  Compressed FASTA is only supported in BGZF format. Use the samtools bgzip utility (instead of gzip) to compress your FASTA.
  ---
  --- Amplicon #2
  Target name	IRX5
  Compressed FASTA is only supported in BGZF format. Use the samtools bgzip utility (instead of gzip) to compress your FASTA.
  ---
  ```
- Get bgzip utility from samtools [htslib-1.10.2](https://github.com/samtools/htslib/releases/download/1.10.2/htslib-1.10.2.tar.bz2) and install
  ```
  cd htslib-1.10.2
  ./configure
  sudo make install
  ```
- Compress file using bgzip
  ```
  bgzip < /Volumes/CRUKCI_HDCLONE/Homo_sapiens.GRCh38.dna.toplevel.fa > Homo_sapiens.GRCh38.dna.toplevel.fa.gz
  ```
- Re-run step 4 with new compressed file: the first time, [pyfaidx](https://pypi.org/project/pyfaidx/) will generate indexed files .fai and it will take some time! There is also a `read_ahead` option to reduce runtime that could be explored.
  ```
  geneditid_create_amplicount_config --project=GEP00009 --genome=/Users/pajanne/workspace/GenEditID/data/reference/Homo_sapiens.GRCh38.dna.toplevel.fa.gz
  ```
- Install [Git Large File Storage (LFS)](https://git-lfs.github.com/) to push reference files to github
  ```
  brew install git-lfs
  git lfs install
  git lfs track "*.fa.gz"
  git lfs track "*.fa.gz.fai"
  git add .gitattributes
  git add data/reference/.
  git commit -m '.....'
  git push
  ```

# Dash and Flask

- [Flask](https://flask.palletsprojects.com/en/1.1.x/)
- [Dash](https://dash.plotly.com/)
- [How to embed a Dash app into an existing Flask app](https://medium.com/@olegkomarov_77860/how-to-embed-a-dash-app-into-an-existing-flask-app-ea05d7a2210b)
