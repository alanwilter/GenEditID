import sys
import os
from operator import itemgetter
from Bio import Align

if len(sys.argv) <= 1:
    print('Please, provide one read counts file (output of count_amplicons.py tool).')
    sys.exit()

with open(sys.argv[1]) as data:
    read_data = False
    results = []
    for line in data:
        if line.startswith('WT amplicon with primers'):
            values = line.split()
            wt_seq = values[-1]
        if read_data:
            line_data = line.strip().split(', ')
            results.append({header[0]: line_data[0],
                            header[1]: line_data[1],
                            header[2]: line_data[2],
                            header[3]: line_data[3],
                            header[4]: line_data[4],
                            header[5]: line_data[5],
                            header[6]: line_data[6],
                            header[7]: line_data[7]})
        if line.startswith('slxid'):
            read_data = True  # start reading data
            header = line.strip().split(', ')

aligner = Align.PairwiseAligner()
aligner.mode = 'local'
sorted_results = sorted(results, key=itemgetter('barcode'))

dirname, basename = os.path.split(sys.argv[1])
align_filename = basename.split('.')[0] + '_align.txt'
fasta_filename = basename.split('.')[0] + '_seq.fa'
guide_loc = 108

with open(align_filename, "w") as out_align, open(fasta_filename, "w") as out_fasta:
    seq_per_barcode = 0
    barcode = ''
    for r in sorted_results:
        if r['barcode'] == barcode:
            seq_per_barcode += 1
        else:
            seq_per_barcode = 1
        seq_header = ">{}_{}_S{}_{}_{}".format(r['slxid'], r['barcode'], seq_per_barcode, r['seq_type'], r['read_count'])
        out_align.write("{}\n".format(seq_header))
        out_fasta.write("{}\n".format(seq_header))
        out_fasta.write("{}\n".format(r['sequence']))
        alignments = aligner.align(wt_seq, r['sequence'])
        #out_align.write("{}\n".format(alignments[0]))
        three_lines_alignment = str(alignments[0]).split()
        print(seq_header)
        print(alignments[0])
        out_align.write("{}\n".format(three_lines_alignment[0][:guide_loc-20].lower() + three_lines_alignment[0][guide_loc-20:guide_loc+20].upper()+three_lines_alignment[0][guide_loc+20:].lower()))
        out_align.write("{}\n".format(three_lines_alignment[1]))
        out_align.write("{}\n".format(three_lines_alignment[2][:guide_loc-20].lower() + three_lines_alignment[2][guide_loc-20:guide_loc+20].upper()+three_lines_alignment[2][guide_loc+20:].lower()))
