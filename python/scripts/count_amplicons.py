import os
import gzip
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import SeqIO
from Bio import Entrez

# ------------------------------------------------------------------------------
# Input data
# search gene id on https://www.ncbi.nlm.nih.gov/nuccore/
# ------------------------------------------------------------------------------

# input data for GEP00004
# submitted_gene_id = 'NC_000006.12'  # ref for chromosom 6!
# submitted_fprimer_seq = 'CTTTTATGTCTTCACAGGGCAGC'  # fprimer_seq
# submitted_rprimer_seq = 'GGTACTGCATTCCGGCCAT'  # = rprimer_seq, reverse complement
#gene_id = 'NG_032093'
#rprimer_seq = Seq('GGTACTGCATTCCGGCCAT', IUPAC.unambiguous_dna)
#rprimer_seq = str(rprimer_seq.reverse_complement())
#fprimer_seq = 'CTTTTATGTCTTCACAGGGCAGC'
#fastq_dir = '/Users/pajon01/mnt/scratchb/genome-editing/GEP00004/fastq/'

# input data for GEP00006
# submitted_gene_id = 'NC_000020.11'  # ref for chromosom 20!
# submitted_fprimer_seq = 'CTTGAGGTTGAAGGTGTTGC'  # = rpimer_seq, reverse complement
# submitted_rprimer_seq = 'GAGCTGCAGAACGGGAACT'  # = fprimer_seq
#gene_id = 'NG_007496'
#rprimer_seq = 'GCAACACCTTCAACCTCAAG'
#fprimer_seq = 'GAGCTGCAGAACGGGAACT'
#fastq_dir = '/Users/pajon01/mnt/scratchb/genome-editing/GEP00006/fastq/'

# input data for GEP00007
# submitted_gene_id = 'ARID1A'
# submitted_fprimer_seq = 'GCAAGATGAGACCTCAGCCA'  # fprimer_seq
# submitted_rprimer_seq = 'TGGACAGCAACAAGGGTCAA'  # = rprimer_seq, reverse complement
gene_id = 'NG_029965'
rprimer_seq = Seq('TGGACAGCAACAAGGGTCAA', IUPAC.unambiguous_dna)
rprimer_seq = str(rprimer_seq.reverse_complement())
fprimer_seq = 'GCAAGATGAGACCTCAGCCA'
fastq_dir = '/Users/pajon01/mnt/scratchb/genome-editing/GEP00007/fastq/'

fastq_extension = '.r_1.fq.gz'  # logic need to change for r_2
frequency_threshold = 5

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
    print('WT amplicon with primers: {} {}'.format(len(wt_amplicon_seq), wt_amplicon_seq))
    print('--------')

# ------------------------------------------------------------------------------
# Occurences of sequences between primers in FASTQ files
# ------------------------------------------------------------------------------
print("slxid, barcode, read_count, read_frequency, seq_type, sequence")
for filename in os.listdir(fastq_dir):
    if filename.endswith(fastq_extension):
        splited_filename = filename.split('.')
        slxid, barcode = splited_filename[0], splited_filename[1]
        occurrences_dict = {}
        reads = 0
        with gzip.open(os.path.join(fastq_dir, filename), "rt") as handle:
            for title, seq, qual in FastqGeneralIterator(handle):
                # for a large FASTQ file, better not use the high-level SeqIO interface
                # for record in SeqIO.parse(handle, "fastq") but work directly with strings
                # to work with quality score converted to integers: record.letter_annotations["phred_quality"]
                fprimer_pos = seq.find(fprimer_seq)
                rprimer_pos = seq.find(rprimer_seq) + len(rprimer_seq)
                amplicon_seq = str(seq)[fprimer_pos:rprimer_pos]
                # find wt sequence
                if amplicon_seq:
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
            print("{}, {}, {}, {:.2f}, {}, {}".format(slxid, barcode, i[0], i[1], seq_type, i[2]))
