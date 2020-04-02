# Things done...

- FastQ files added into GDrive to be [downloaded from here](https://drive.google.com/drive/folders/1MN_vzy3hjGOAnycwwtI53nrAuOaB5RJf?usp=sharing) (step 3)
- Access to [Google Drive paper from here](https://drive.google.com/drive/folders/1MQAmhxjuewH2gDoUkzXzz1wmgMK6CV7E?usp=sharing)
- Slack integration, find us on [geneditid](geneditid.slack.com)

# Things in progress...

- Create proper test data
  - [x] select 3 FastQ files
  - [ ] run all analysis steps on these 3 files only

# Things to do...

- Simplify submission spreadsheet and database
- Add new plots to the WebApp
  - Add output of the run_ampli_count tool into a new database table for plotting (add class AmpliCountResult in model.py)
  - Update the plotting scripts to call the code from the WebApp to avoid code duplication

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

  - Check results in `amplicount.csv`
