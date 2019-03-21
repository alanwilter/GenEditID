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

- DONE. update gene ids of all project to be Ensembl ones
```
python python/scripts/update_target_gene_ids.py
```

- DONE. amplifind: find amplicon and target coordinates using pyfaidx
and modify `create_pipeline_files.py` to use amplifind
and check amplicon and target coordinates of projects GEP00005 / GEP00009 and GEP00010

```
python python/scripts/create_pipeline_files.py --project=GEP00010 --genome=/Users/pajon01/mnt/refdata/reference_genomes/homo_sapiens/GRCh38_hs38d1/fasta/hsa.GRCh38_hs38d1.fa --seq-dict=/Users/pajon01/mnt/refdata/reference_genomes/homo_sapiens/GRCh38_hs38d1/fasta/hsa.GRCh38_hs38d1.dict --filelist=/Users/pajon01/mnt/scratchb/genome-editing/GEP00010/filelist.csv
diff targets.txt /Users/pajon01/mnt/scratchb/genome-editing/GEP00010/targets.txt
diff amplicons.txt /Users/pajon01/mnt/scratchb/genome-editing/GEP00010/amplicons.txt

python python/scripts/create_pipeline_files.py --project=GEP00009 --genome=/Users/pajon01/mnt/refdata/reference_genomes/homo_sapiens/GRCh38_hs38d1/fasta/hsa.GRCh38_hs38d1.fa --seq-dict=/Users/pajon01/mnt/refdata/reference_genomes/homo_sapiens/GRCh38_hs38d1/fasta/hsa.GRCh38_hs38d1.dict --filelist=/Users/pajon01/mnt/scratchb/genome-editing/GEP00009/filelist.csv
diff targets.txt /Users/pajon01/mnt/scratchb/genome-editing/GEP00009/targets.txt
diff amplicons.txt /Users/pajon01/mnt/scratchb/genome-editing/GEP00009/amplicons.txt

python python/scripts/create_pipeline_files.py --project=GEP00005 --genome=/Users/pajon01/mnt/refdata/reference_genomes/homo_sapiens/GRCh38_hs38d1/fasta/hsa.GRCh38_hs38d1.fa --seq-dict=/Users/pajon01/mnt/refdata/reference_genomes/homo_sapiens/GRCh38_hs38d1/fasta/hsa.GRCh38_hs38d1.dict --filelist=/Users/pajon01/mnt/scratchb/genome-editing/GEP00005v1/filelist.csv
diff targets.txt /Users/pajon01/mnt/scratchb/genome-editing/GEP00005v1/targets.txt
diff amplicons.txt /Users/pajon01/mnt/scratchb/genome-editing/GEP00005v1/amplicons.txt

```
- DONE. re-run amplicon pipeline for projects GEP00005 / GEP00009 and GEP00010 with new coordinates
finally solved the issues with the coordinates by using pyfaidx against fasta reference genome and guide location +/- 1000 bases and primer pair

- DONE. amplicount:
  - read distribution per amplicon:
    - table summary: project id . barcode . total reads . amplicon1 . amplicon2 ...
  - read distribution per allele per amplicon:
    - table: barcode . amplicon1 . edit1 . ge or wt . seq
  - re-run read counts on projects GEP00005 / GEP00009 and GEP00010

- ampliplot:
  - DONE. plot summary per amplicon: bar plot for reads per amplicon for all samples
  - plot bar plot per allele for all samples

- ampliscore:
  - read 2 csv output of amplicon pipelines
  - read 1 csv output of amplicount
  - calculate score
  - plot score across all samples

- amplifind: make it work without db only using input csv file
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

- install http://bioconductor.org/packages/release/bioc/html/amplican.html and run on GEP00005 / GEP00009 and GEP00010
- install https://github.com/lucapinello/CRISPResso and run on GEP00005 / GEP00009 and GEP00010

- Exploring variant annotation with VEP
curl 'https://rest.ensembl.org/vep/human/region/9:22125502-22125503:1/C?' -H 'Content-type:application/json'
