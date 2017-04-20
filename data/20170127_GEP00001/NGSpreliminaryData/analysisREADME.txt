All the data on cluster 1
Path : /mnt/scratcha/bioinformatics/chilam01/projects/CRISPR/NGS/20150522_ChanD_DF_CRISPR/SLX-13775

Scripts : chanduAnalysis
fastq : fastq
fastQC output : fastq/fastqc

look at the README.txt for other steps



Why some samples showed low alignament rate?

alignment report 'SLX-13775.alignment_coverage_report.html' showed that alignmnet rate varies widely, ranging from 48.7 % in sample FLD0226 to 99.5% in the sample FLD0194.

I investigated the why some samples showed low alignmnet rate, by taking FLD0226 data in the first instance.

Total number of reads 120296 ( R1 60148 + R2 60148 )

Read length = 150

get unaligned reads from bam file (samtools view -f 4 ../../FLD0226.bwamem.bam )
samtools view -f 4 FLD0226.bwamem.bam  > alnStats/FLD0226_unaligned.sam

Get unaligned read ids
cat  FLD0226_unaligned.sam | cut -f 1 | sort | uniq >FLD0226_unaligned_read_ids

Copy read 1 fastq file
cp /mnt/scratcha/bioinformatics/chilam01/projects/CRISPR/NGS/20150522_ChanD_DF_CRISPR/SLX-13775/fastq/SLX-13775.FLD0226.000000000-B4LYY.s_1.r_1.fq.gz .

gunzip SLX-13775.FLD0226.000000000-B4LYY.s_1.r_1.fq.gz 


Using 'getUnalignedReadsInFastaFormat_v1.0.R' script, I extracted unaligned reads and saved them as both as fasta and fastq files.

fastq file is for fastQC and fasta is for pairwise alignmenmt.

mkdir FLD0226_unaligned_fastqc

run factqc:
 /Users/chilam01/Apps/fastQC/0.11.5/FastQC/fastqc FLD0226_unaligned_reads.fastq -o FLD0226_unaligned_fastqc --extract
 
 factQC report showed significant drop in base qualities after 90 bp.
 
 
 
 # pairwise alignment
 
 FLD0226 - amplicon - STAT3.2.in__STAT3.3.in__STAT3.4.in
 
 Get gequence from this regions.
 samtools faidx /mnt/scratchb/bioinformatics/reference_data/reference_genomes/homo_sapiens/GRCh37_g1kp2/fasta/hsa.GRCh37_g1kp2.fa chr17:40497540-40497740 >STAT3.2.in__STAT3.3.in__STAT3.4.in_ref_seq.txt 
 
 
 run pairwise alignment by using exonerate
 
 sabtch run_pairwiseAln_v2.0.sh
 
 
 
 There are total 32,569 unaligned reads, but more than 25,000 reads aligned only first 80 bases with indels.
 
 Questions?
 
 - Why last 70 based not aligning, is it beacuse of bad quality or something else?
 - About 4,000 reads showed popy T, ( *.T.*), whay is that?
 - After 90 there is a poly A sequence in most of the sequences (look at FLD0226_unaligned_only_seq), why is that?
 
 
 # Dose read2 show same patteren?
 
 cp ../../../fastq/SLX-13775.FLD0226.000000000-B4LYY.s_1.r_2.fq.gz .
 gunzip SLX-13775.FLD0226.000000000-B4LYY.s_1.r_2.fq.gz 
 
source( file='getUnalignedReadsInFastaFormat_v2.0.R', echo=T)

sbatch run_pairwiseAln_v3.0.sh 


# fastQC

mkdir FLD0226_unaligned_fastqc_r2


 /Users/chilam01/Apps/fastQC/0.11.5/FastQC/fastqc FLD0226_unaligned_reads_r2.fastq -o FLD0226_unaligned_fastqc_r2 --extract
 
 ~/bin/MakeFastaToSingleLineFasta.pl FLD0226_unaligned_reads_r2.fasta  | grep -v '^>' >FLD0226_unaligned_only_seq_r2
 


# read1 and read 2 showed similar patters


# FLD0194 is  good sample showed 99.5% alignment rate and 99.5% useble bases and 97% asigned read pairs.

How this data looks like, when compared to FLD0226?
cp ../../../fastq/SLX-13775.FLD0194.000000000-B4LYY.s_1.r_* .
gunzip *

# fastq to fasta
/home/chilam01/Apps/fastx_toolkit/fastx_toolkit-0.0.14/bin/fastq_to_fasta -i SLX-13775.FLD0194.000000000-B4LYY.s_1.r_1.fq -o SLX-13775.FLD0194.000000000-B4LYY.s_1.r_1.fasta

/home/chilam01/Apps/fastx_toolkit/fastx_toolkit-0.0.14/bin/fastq_to_fasta -i SLX-13775.FLD0194.000000000-B4LYY.s_1.r_2.fq -o SLX-13775.FLD0194.000000000-B4LYY.s_1.r_2.fasta

cat SLX-13775.FLD0194.000000000-B4LYY.s_1.r_1.fasta | grep -v '^>' > SLX-13775.FLD0194.000000000-B4LYY.s_1.r_1_only_seq

cat SLX-13775.FLD0194.000000000-B4LYY.s_1.r_2.fasta | grep -v '^>' > SLX-13775.FLD0194.000000000-B4LYY.s_1.r_2_only_seq
