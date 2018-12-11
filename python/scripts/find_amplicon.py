from Bio.Seq import Seq
from Bio import SeqIO
from Bio import Entrez

# input data for GEP00012
gene_id = 'NM_017628'
gene_start = 105145875
rprimer_seq = 'GGAGACACCAAGTGGCACTC'
fprimer_seq = 'TCAACTAGAGGGCAGCCTTG'

Entrez.email = "A.N.Other@example.com"
handle = Entrez.esearch(db="nucleotide", term=gene_id)
record = Entrez.read(handle)
handle = Entrez.efetch(db="nucleotide", id=record["IdList"], rettype="gb", retmode="text")
for record in SeqIO.parse(handle, "gb"):
    print('--------')
    print('Gene ID: {}'.format(record.name))
    print('Forward primer: {} {}'.format(len(fprimer_seq), fprimer_seq))
    print('Reverse primer: {} {}'.format(len(rprimer_seq), rprimer_seq))
    fprimer_start_pos = record.seq.find(fprimer_seq)
    rprimer_start_pos = record.seq.find(rprimer_seq)
    target_start = gene_start
    target_end = gene_start + rprimer_start_pos + len(rprimer_seq)
    print('Target coordinates: {} {}'.format(target_start, target_end))
    amplicon_start = gene_start + fprimer_start_pos + len(fprimer_seq)
    amplicon_end = gene_start + rprimer_start_pos
    print('Amplicon coordinates: {} {}'.format(amplicon_start, amplicon_end))
    print('--------')
