# kickstart
sbatch runKickStart_v1.0.sh 

mkdir bwa

cd bwa
cp /mnt/scratcha/bioinformatics/chilam01/projects/CRISPR/NGS/20150522_ChanD_DF_CRISPR/bwa/create_bwa_config.sh .

 # edit reference genome paths

# create alignment config file
sh create_bwa_config.sh SLX-13775 /mnt/scratcha/bioinformatics/chilam01/projects/CRISPR/NGS/20150522_ChanD_DF_CRISPR/SLX-13775/fastq general alignmentConfig_v1.0.xml

# queue changed from bioinformatics to general in alignmentConfig_v1.0.xml

cp /mnt/scratcha/bioinformatics/chilam01/projects/CRISPR/NGS/20150522_ChanD_DF_CRISPR/bwa/runAlignment_v1.0.sh .

# edit runAlignment_v1.0.sh 
# add this line

SOFTWARE_ROOT=/home/bioinformatics/pipelinesoftware/alignment/el7
sbatch runAlignment_v1.0.sh



# prepare amplicon and target files
cd /mnt/scratcha/bioinformatics/chilam01/projects/CRISPR/NGS/20150522_ChanD_DF_CRISPR/SLX-13775/refDataDir

Rscript createAmpliconAndTargetFiles_v3.0.R

cat matt/sequence_dictionary.txt amplicons_v1.txt >amplicons_v2.txt 
cat matt/sequence_dictionary.txt target_v1.txt >target_v2.txt 

pwd : /mnt/scratcha/bioinformatics/chilam01/projects/CRISPR/NGS/20150522_ChanD_DF_CRISPR/SLX-13775

cp /home/bioinformatics/software/pipelines/ampliconseq/bin/create_samples_file.sh .

sh create_samples_file.sh samplesheet.csv samples.txt

cp /home/bioinformatics/software/pipelines/ampliconseq/bin/create_ampliconseq_config.sh .

# Above script did not work for me

# copy config file from Matt run folder and modify it
cp /mnt/scratcha/bioinformatics/eldrid01/20170131_PiskorzA_JB_ampliconseq/SLX-10379/config.xml .

# run amplicon seq pipeline
#I used perlbrew to intall perl locally and cpam to install bioperl packages
# before starting the pipeline I need specify to use local perl
perlbrew use perl-5.16.0
sbatch runAmpliconseq_vardict_v1.0.sh 

# haplotype caller
sbatch runAmpliconseq_HaplotypeCaller_v1.0.sh 

