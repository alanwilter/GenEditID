# List of things to do

- DONE. check fastq-join quality on project GEP00013 (automatised the task and gather output)

```
ssh clust1-headnode
cd /scratchb/bioinformatics/pajon01/genome-editing/GEP00013/fastq
for f in *.s_1.r_1.fq.gz; do echo "~/fastq-join/fastq-join $f `echo $f | cut -d'.' -f1-4`.r_2.fq.gz -o `echo $f | cut -d'.' -f1-4`.fq"; done

~/fastq-join/fastq-join -h
Usage: fastq-join [options] <read1.fq> <read2.fq> [mate.fq] -o <read.0q>
Version: 1.3.1

Joins two paired-end reads on the overlapping ends.

Options:

-o FIL     See 'Output' below
-v C       Verifies that the 2 files probe id's match up to char C
            use ' ' (space) for Illumina reads
-p N       N-percent maximum difference (8)
-m N       N-minimum overlap (6)
-r FIL     Verbose stitch length report
-R         No reverse complement
-x         Allow insert < read length
```

```
cp shell/ngs/job_joinreads.sh /Users/pajon01/mnt/scratchb/genome-editing/GEP00013/.
```

```
ssh clust1-headnode
cd /scratchb/bioinformatics/pajon01/genome-editing/GEP00013/
sbatch job_joinreads.sh
cat joinreads.out
```
- DONE. create similar script for merging reads (automatised the task and gather output)

```
scp shell/ngs/job_mergereads.sh clust1-headnode:/scratchb/bioinformatics/pajon01/genome-editing/GEP00010/.
ssh clust1-headnode
cd /scratchb/bioinformatics/pajon01/genome-editing/GEP00010
sbatch job_mergereads.sh
tail -f mergereads.out
```

- DONE. get 2bit file for GRCh38_hs38d1 using faToTwoBit in.fa out.2bit

```
ssh clust1-headnode
cd /home/pajon01/bin
ln -s /home/bioinformatics/software/ucsc-tools/ucsc-tools-20171107/faToTwoBit
cd /scratchb/bioinformatics/pajon01/genome-editing/check_amplicon_coords
ll /scratcha/bioinformatics/reference_data/reference_genomes/homo_sapiens/GRCh38_hs38d1/fasta/hsa.GRCh38_hs38d1.fa
sbatch job_faToTwoBit.sh
```

- DONE. update gene ids of all project to be Ensembl ones
```
python python/scripts/update_target_gene_ids.py
```

- DONE. run_ampli_count:
  - read distribution per amplicon:
    - table summary: project id . barcode . total reads . amplicon1 . amplicon2 ...
  - read distribution per allele per amplicon:
    - table: barcode . amplicon1 . edit1 . ge or wt . seq
  - re-run read counts on projects GEP00005 / GEP00009 and GEP00010

- run_variant_id:
  - DONE. plot summary per amplicon: bar plot for reads per amplicon for all samples
  - DONE. plot bar plot per allele for all samples

- run_variant_id:
  - DONE. calculate score
  - DONE. plot score across all samples

- run_ampli_find: make it work without db only using input csv file
  - refgenome fasta file
  - config file: guide_loc,chr,strand,fprimer,rprimer
- add read count table into db and webapp
- simplify submission form
- primerfind: integrate blat to find primer pairs in genome
```
./gfServer start localhost 8888 GRCh38_hs38d1.2bit -log=gfServer.log -canStop -stepSize=5 > gfServer.out &
./gfServer status localhost 8888
./gfServer pcr localhost 8888 GGGCCCCAGAAATTCAAAGG AGAGCAGGAGGTCCAGTCAG 300
GRCh38_hs38d1.2bit:chr12	57600720	57600945	+
```
