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

- DONE. run amplicon pipeline on GEP00013
- DONE. load incucyte data files for GEP00013
- load incucyte 1in10 data files for GEP00013 (add plates)

- update gene ids of all project to be Ensembl ones
```
GEP00012 TET2 chr4 NM_017628 ENSG00000168769 ERR: 2 identical targets/primers
GEP00011 ETV6 chr12 NM_001987.4 ENSG00000139083 ERR: 2 identical targets/primers
GEP00010 FTO NC_000016.10 ENSG00000140718
GEP00010 LEPR NC_000001.11 ENSG00000116678
GEP00010 CCKAR NC_000004.12 ENSG00000163394
GEP00009 FTO NC_000016.10 ENSG00000140718
GEP00009 IRX5 NC_000016.10 ENSG00000176842
GEP00007 ARID1A NC_000001.11 ENSG00000117713 ERR: 2 identical targets/primers
GEP00006 Jag1 NC_000020.11 ENSG00000101384 ERR: 2 identical targets/primers

GEP00005 FTO NC_000016.10 ENSG00000140718
GEP00005 RPGRIP1L NC_000016.10 ENSG00000103494
GEP00005 PREPL NC_000002.12 ENSG00000138078
GEP00005 IRX3 NC_000016.10 ENSG00000177508
GEP00005 IRX5 NC_000016.10 ENSG00000176842
GEP00005 LEPRY2 NC_000001.11 ENSG00000116678 single base modification
GEP00005 LEPRY3 NC_000001.11 ENSG00000116678 single base modification
GEP00005 LEPRY4 NC_000001.11 ENSG00000116678 single base modification
GEP00005 LEPRLF NC_000001.11 ENSG00000116678 single base modification

GEP00004 ARID1B NC_000006.12 ENSG00000049618 ERR: 2 identical targets/primers
GEP00003 GFP GFP ???
GEP00002 PTEN chr10 NC_000010.10 ENSG00000171862
GEP00002 PTEN NC_000010.10
GEP00002 PTEN NC_000010.10
GEP00002 PTEN NC_000010.10

```

- integrate blat to find amplicon and target coordinates (use GRCh38_hs38d1 instead of GRCh38)

```
./gfServer start localhost 8888 GRCh38_hs38d1.2bit -log=gfServer.log -canStop -stepSize=5 > gfServer.out &
./gfServer status localhost 8888
./gfServer pcr localhost 8888 GGGCCCCAGAAATTCAAAGG AGAGCAGGAGGTCCAGTCAG 300
GRCh38_hs38d1.2bit:chr12	57600720	57600945	+
```

or integrate pyfaidx using amplicon start +/-500

- find amplicon and target coordinates of projects GEP00005 / GEP00009 and GEP00010
- automatise read counts to run on multiple amplicons (run projects GEP00005 / GEP00009 and GEP00010)
- re-run amplicon pipeline for projects GEP00005 / GEP00009 and GEP00010
- install https://github.com/lucapinello/CRISPResso and run on GEP00005 / GEP00009 and GEP00010
- install http://bioconductor.org/packages/release/bioc/html/amplican.html and run on GEP00005 / GEP00009 and GEP00010
