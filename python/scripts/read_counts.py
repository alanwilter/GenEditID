import os
import sys
import yaml
import gzip
import argparse
import log as logger
from Bio import SeqIO
from dnascissors.utils import get_target_sequence


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("--config", dest="config", action="store", help="The yaml config file e.g. 'amplicons.yml'", required=True)
    parser.add_argument("--output", dest="out", action="store", help="The output file e.g. 'GEPID_read_counts.csv'.", required=True)
    options = parser.parse_args()

    log = logger.get_custom_logger(os.path.join(os.path.dirname(__file__), 'read_counts.log'))

    with open(options.config, 'r') as yml_file:
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
    # Retrieve target sequence from Ensembl
    # ------------------------------------------------------------------------------
    fprimer_pos, fprimer_seq, rprimer_pos, rprimer_seq, wt_amplicon_seq = get_target_sequence(ensembl_gene_id=gene_id, fprimer_seq=fprimer_seq, rprimer_seq=rprimer_seq)

    header = ('--------\n')
    header += ('Gene id: {}\n'.format(gene_id))
    header += ('Forward primer: {} {}\n'.format(len(fprimer_seq), fprimer_seq))
    header += ('Reverse primer: {} {}\n'.format(len(rprimer_seq), rprimer_seq))
    if fprimer_pos == -1 or rprimer_pos == -1:
        header += ('Primers not found in Ensembl Gene ID {}!\n'.format(ensembl_gene_id))
        log.warning(header)
        sys.exit()
    else:
        header += ('WT amplicon sequence: {} {}\n'.format(len(wt_amplicon_seq), wt_amplicon_seq))
        header += ('--------')

    # ------------------------------------------------------------------------------
    # Occurences of sequences between primers in FASTQ files
    # ------------------------------------------------------------------------------
    header += ('--------\n')
    header += ('FastQ files: {}\n'.format(fastq_dir))
    header += ('FastQ extension: {}\n'.format(fastq_extension))
    header += ('Read quality Threshold: {}\n'.format(quality_threshold))
    header += ('Frequency Threshold: {} %\n'.format(frequency_threshold))
    header += ('Count reads with a minimum base quality score above the read quality threshold within the amplicon.\n')
    header += ('Report only sequences with read counts which are above the frequency threshold for this amplicon.\n')
    header += ('--------\n')
    header += ("slxid, barcode, total_reads, total_filtered_reads, read_count, read_frequency, seq_type, sequence\n")
    log.info(header)

    with open(options.out, 'w') as out:
        out.write(header)
        for filename in os.listdir(fastq_dir):
            if filename.endswith(fastq_extension):
                splited_filename = filename.split('.')
                slxid, barcode = splited_filename[0], splited_filename[1]
                occurrences_dict = {}
                total_reads = 0
                total_filtered_reads = 0
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
                            if (fprimer_pos > rprimer_pos):
                                fprimer_pos = 0
                            if min(rec.letter_annotations["phred_quality"][fprimer_pos:rprimer_pos]) >= quality_threshold:
                                occurrences_dict[amplicon_seq] = occurrences_dict.get(amplicon_seq, 0) + 1
                                total_filtered_reads += 1

                occurrences_items = [(((occ / total_filtered_reads) * 100), occ, seq) for seq, occ in occurrences_dict.items() if ((occ / total_filtered_reads) * 100) > frequency_threshold]
                occurrences_items.sort()
                occurrences_items.reverse()
                occurrences_items = [(occ, f, seq) for f, occ, seq in occurrences_items]
                for i in occurrences_items:
                    seq_type = 'ge'
                    if i[2] == (fprimer_seq + wt_amplicon_seq + rprimer_seq):
                        seq_type = 'wt'
                    out.write("{}, {}, {}, {}, {}, {:.2f}, {}, {}\n".format(slxid, barcode, total_reads, total_filtered_reads, i[0], i[1], seq_type, i[2]))


if __name__ == '__main__':
    main()
