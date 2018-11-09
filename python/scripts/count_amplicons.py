import os
import yaml
import gzip
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import SeqIO
from Bio import Entrez

yml_filepath = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'amplicons.yml')

with open(yml_filepath, 'r') as yml_file:
    amplicons = yaml.load(yml_file)

# ------------------------------------------------------------------------------
# Input data in config file amplicons.yml
# ------------------------------------------------------------------------------
gene_id = amplicons['gene_id']
rprimer_seq = amplicons['rprimer_seq']
rprimer_seq = str(Seq(amplicons['rprimer_seq'], IUPAC.unambiguous_dna).reverse_complement()) # for all other project than 12
#rprimer_seq = amplicons['rprimer_seq'] # for project 12
fprimer_seq = amplicons['fprimer_seq']
fastq_dir = amplicons['fastq_dir']

fastq_extension = amplicons['fastq_extension']
quality_threshold = amplicons['quality_threshold']
frequency_threshold = amplicons['frequency_threshold']

# ------------------------------------------------------------------------------
# Retrieve wild type amplicon sequence from Nucleotide database on NCBI
# ------------------------------------------------------------------------------
Entrez.email = "A.N.Other@example.com"
handle = Entrez.esearch(db="nucleotide", term=gene_id)
record = Entrez.read(handle)
handle = Entrez.efetch(db="nucleotide", id=record["IdList"], rettype="gb", retmode="text")
for record in SeqIO.parse(handle, "gb"):
    print('--------')
    print('Gene ID: {}'.format(record.name))
    print('Forward primer: {} {}'.format(len(fprimer_seq), fprimer_seq))
    print('Reverse primer: {} {}'.format(len(rprimer_seq), rprimer_seq))
    fprimer_pos = record.seq.find(fprimer_seq)
    rprimer_pos = record.seq.find(rprimer_seq) + len(rprimer_seq)
    wt_amplicon_seq = str(record.seq)[fprimer_pos:rprimer_pos]
    #wt_amplicon_seq = wt_amplicon_seq_whole[:96] + wt_amplicon_seq_whole[135:]
    print('WT amplicon with primers: {} {}'.format(len(wt_amplicon_seq), wt_amplicon_seq))
    print('--------')

# ------------------------------------------------------------------------------
# Occurences of sequences between primers in FASTQ files
# ------------------------------------------------------------------------------
print('--------')
print('FastQ files: {}'.format(fastq_dir))
print('Read quality Threshold: {}'.format(quality_threshold))
print('Frequency Threshold: {} %'.format(frequency_threshold))
print('Counting reads which have both primers in their sequence, and with a minimum base quality score above the threshold within the amplicon sequence.')
print('Reporting only sequences with read counts which are above the threshold.')
print('--------')
print("slxid, barcode, total_reads, total_filtered_reads, read_count, read_frequency, seq_type, sequence")
for filename in os.listdir(fastq_dir):
    if filename.endswith(fastq_extension):
        splited_filename = filename.split('.')
        slxid, barcode = splited_filename[0], splited_filename[1]
        occurrences_dict = {}
        total_reads = 0
        reads = 0
        # filtering reads based on quality
        with gzip.open(os.path.join(fastq_dir, filename), "rt") as handle:
            for rec in SeqIO.parse(handle, "fastq"):
                total_reads += 1
                fprimer_pos = str(rec.seq).find(fprimer_seq)
                rprimer_pos = str(rec.seq).find(rprimer_seq)
                # check both primers found
                if not fprimer_pos == -1 and not rprimer_pos == -1:
                    rprimer_pos = rprimer_pos + len(rprimer_seq)
                    amplicon_seq = str(rec.seq)[fprimer_pos:rprimer_pos]
                    #amplicon_seq = amplicon_seq_whole[:96] + amplicon_seq_whole[135:]
                    # check min base quality is above 10
                    # check quality over a window instead of full sequence
                    if min(rec.letter_annotations["phred_quality"][fprimer_pos:rprimer_pos]) >= quality_threshold:
                        occurrences_dict[amplicon_seq] = occurrences_dict.get(amplicon_seq, 0) + 1
                        reads += 1

        occurrences_items = [(((occ / reads) * 100), occ, seq) for seq, occ in occurrences_dict.items() if ((occ / reads) * 100) > frequency_threshold]
        occurrences_items.sort()
        occurrences_items.reverse()
        occurrences_items = [(occ, f, seq) for f, occ, seq in occurrences_items]
        for i in occurrences_items:
            seq_type = 'ge'
            if i[2] == wt_amplicon_seq:
                seq_type = 'wt'
            print("{}, {}, {}, {}, {}, {:.2f}, {}, {}".format(slxid, barcode, total_reads, reads, i[0], i[1], seq_type, i[2]))
