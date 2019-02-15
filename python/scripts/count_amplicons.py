import os
import sys
import yaml
import gzip
import requests
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC


yml_filepath = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'amplicons.yml')

with open(yml_filepath, 'r') as yml_file:
    amplicons = yaml.load(yml_file)

# ------------------------------------------------------------------------------
# Input data in config file amplicons.yml
# ------------------------------------------------------------------------------
gene_id = amplicons['gene_id']
fprimer_seq = amplicons['fprimer_seq']
rprimer_seq = amplicons['rprimer_seq']

fastq_dir = amplicons['fastq_dir']
fastq_extension = amplicons['fastq_extension']
quality_threshold = amplicons['quality_threshold']
frequency_threshold = amplicons['frequency_threshold']

# ------------------------------------------------------------------------------
# Retrieve wild type amplicon sequence from Ensembl
# ------------------------------------------------------------------------------
url = ("http://rest.ensembl.org/sequence/id/{}?").format(gene_id)
r = requests.get(url, headers={"Content-Type": "text/plain"})
sequence = r.text
fprimer_pos = sequence.find(fprimer_seq)
rprimer_seq = str(Seq(rprimer_seq, IUPAC.unambiguous_dna).reverse_complement())
rprimer_pos = sequence.find(rprimer_seq)
wt_amplicon_seq = sequence[fprimer_pos+len(fprimer_seq):rprimer_pos-len(rprimer_seq)]

print('--------')
print('Gene id: {}'.format(gene_id))
print('Forward primer: {} {}'.format(len(fprimer_seq), fprimer_seq))
print('Reverse primer: {} {}'.format(len(rprimer_seq), rprimer_seq))
print('WT amplicon sequence: {} {}'.format(len(wt_amplicon_seq), wt_amplicon_seq))
print('--------')

if fprimer_pos == -1 or rprimer_pos == -1:
    print('Primers not found! [forward={}, reverse={}]'.format(fprimer_pos, rprimer_pos))
    sys.exit()

# ------------------------------------------------------------------------------
# Occurences of sequences between primers in FASTQ files
# ------------------------------------------------------------------------------
print('--------')
print('FastQ files: {}'.format(fastq_dir))
print('FastQ extension: {}'.format(fastq_extension))
print('Read quality Threshold: {}'.format(quality_threshold))
print('Frequency Threshold: {} %'.format(frequency_threshold))
print('Count reads with a minimum base quality score above the read quality threshold within the amplicon.')
print('Report only sequences with read counts which are above the frequency threshold for this amplicon.')
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
                    # check min base quality is above quality_threshold
                    # check quality over a window instead of full sequence
                    if (fprimer_pos > rprimer_pos) and not rprimer_pos == -1:
                        fprimer_pos = 0
                    if min(rec.letter_annotations["phred_quality"][fprimer_pos:rprimer_pos]) >= quality_threshold:
                        occurrences_dict[amplicon_seq] = occurrences_dict.get(amplicon_seq, 0) + 1
                        reads += 1

        occurrences_items = [(((occ / reads) * 100), occ, seq) for seq, occ in occurrences_dict.items() if ((occ / reads) * 100) > frequency_threshold]
        occurrences_items.sort()
        occurrences_items.reverse()
        occurrences_items = [(occ, f, seq) for f, occ, seq in occurrences_items]
        for i in occurrences_items:
            seq_type = 'ge'
            if i[2] == (fprimer_seq + wt_amplicon_seq + rprimer_seq):
                seq_type = 'wt'
            print("{}, {}, {}, {}, {}, {:.2f}, {}, {}".format(slxid, barcode, total_reads, reads, i[0], i[1], seq_type, i[2]))
