import os
import csv
import gzip
import argparse
import collections
import log as logger
from Bio import SeqIO


def count_reads(log, outputfile, fastq_dir, fastq_extension, amplicons, quality_threshold, reads_threshold, frequency_threshold):
    filename, ext = os.path.splitext(outputfile)
    with open(outputfile, 'w') as out, open("{}_coverage{}".format(filename, ext), 'w') as out_coverage:
        out.write("slxid,barcode,total_reads,amplicon_id,amplicon_filtered_reads,variant_reads,variant_frequency,seq_type,sequence\n")
        out_coverage_header = True
        for filename in sorted(os.listdir(fastq_dir)):
            if filename.endswith(fastq_extension):
                splited_filename = filename.split('.')
                slxid, barcode = splited_filename[0], splited_filename[1]
                total_reads = 0
                amplicon_filtered_reads = collections.OrderedDict()
                variant_occurences = collections.OrderedDict()
                for amplicon in amplicons:
                    variant_occurences[amplicon['id']] = {}
                    amplicon_filtered_reads[amplicon['id']] = 0
                if out_coverage_header:
                    out_coverage.write("slxid,barcode,total_reads,{}\n".format(','.join(amplicon_filtered_reads.keys())))
                    out_coverage_header = False
                # filtering reads based on quality
                with gzip.open(os.path.join(fastq_dir, filename), "rt") as handle:
                    log.info('Counting reads in file {}'.format(filename))
                    for rec in SeqIO.parse(handle, "fastq"):
                        total_reads += 1
                        for amplicon in amplicons:
                            fprimer_pos = str(rec.seq).find(amplicon['fprimer'])
                            rprimer_pos = str(rec.seq).find(amplicon['rprimer'])
                            # check both primers found
                            if not fprimer_pos == -1 and not rprimer_pos == -1:
                                rprimer_pos = rprimer_pos + len(amplicon['rprimer'])
                                variant_seq = str(rec.seq)[fprimer_pos:rprimer_pos]
                                # check min base quality is above quality_threshold
                                # check quality over a window instead of full sequence
                                if (fprimer_pos > rprimer_pos):
                                    fprimer_pos = 0
                                if min(rec.letter_annotations["phred_quality"][fprimer_pos:rprimer_pos]) >= int(quality_threshold):
                                    variant_occurences[amplicon['id']][variant_seq] = variant_occurences[amplicon['id']].get(variant_seq, 0) + 1
                                    amplicon_filtered_reads[amplicon['id']] = amplicon_filtered_reads.get(amplicon['id'], 0) + 1
                    log.info("{},{},{},{}".format(slxid, barcode, total_reads, ','.join(str(n) for n in amplicon_filtered_reads.values())))
                    out_coverage.write("{},{},{},{}\n".format(slxid, barcode, total_reads, ','.join(str(n) for n in amplicon_filtered_reads.values())))

                for amplicon_id in variant_occurences:
                    for seq in variant_occurences[amplicon_id]:
                        seq_type = 'var'
                        for amplicon in amplicons:
                            if amplicon_id == amplicon['id'] and seq == amplicon['amplicon']:
                                seq_type = 'REF'
                        variant_reads = variant_occurences[amplicon_id][seq]
                        variant_frequency = (variant_reads / amplicon_filtered_reads[amplicon_id]) * 100
                        if (amplicon_filtered_reads[amplicon_id] >= int(reads_threshold)) and (variant_frequency > int(frequency_threshold)):
                            out.write("{},{},{},{},{},{},{:.2f},{},{}\n".format(slxid, barcode, total_reads, amplicon_id, amplicon_filtered_reads[amplicon_id], variant_reads, variant_frequency, seq_type, seq))


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("--config", dest="config", action="store", help="The 4 columns input config file: 'id,fprimer,rprimer,amplicon'", default='amplicount_config.csv', required=False)
    parser.add_argument("--fastqdir", dest="fastq_dir", action="store", help="Fastq file directory", required=True)
    parser.add_argument("--fastqext", dest="fastq_extension", action="store", help="Fastq file extension", default='.fqjoin.gz', required=False)
    parser.add_argument("--quality", dest="quality_threshold", action="store", help="Threshold for min phred quality across the amplicon sequence", default=10, required=False)
    parser.add_argument("--reads", dest="reads_threshold", action="store", help="Threshold for min number of filtered reads to report variant", default=1000, required=False)
    parser.add_argument("--frequency", dest="frequency_threshold", action="store", help="Frequency threshold to report variant over filtered reads", default=5, required=False)
    parser.add_argument("--output", dest="output", action="store", help="The output file", default='amplicount.csv', required=False)
    options = parser.parse_args()

    log = logger.get_custom_logger('amplicount.log')
    log.info('>>> Getting list of amplicons...')
    amplicons = []
    with open(options.config) as config:
        reader = csv.DictReader(config)
        for row in reader:
            amplicons.append(row)
    for amplicon in amplicons:
        log.info('Amplicon: {}'.format(amplicon['id']))
        log.info('Forward primer: {} {}'.format(len(amplicon['fprimer']), amplicon['fprimer']))
        log.info('Reverse primer: {} {}'.format(len(amplicon['rprimer']), amplicon['rprimer']))
        log.info('Ref amplicon seq: {} {}'.format(len(amplicon['amplicon']), amplicon['amplicon']))
        log.info('----------')

    log.info('>>> Counting reads per variant per amplicon...')
    log.info('FastQ file directory: {}'.format(options.fastq_dir))
    log.info('FastQ extension: {}'.format(options.fastq_extension))
    log.info('Read quality Threshold: {}'.format(options.quality_threshold))
    log.info('Read count Threshold: {}'.format(options.reads_threshold))
    log.info('Frequency Threshold: {} %'.format(options.frequency_threshold))
    log.info('Count reads with a minimum base quality score above the read quality threshold within the amplicon.')
    log.info('Report only variants with amplicon read counts which are above the read count threshold for this amplicon.')
    log.info('and report only variants with variant read counts which are above the frequency threshold for this amplicon.')

    # count reads associated with each variant for each amplicon
    count_reads(log,
                options.output,
                options.fastq_dir,
                options.fastq_extension,
                amplicons,
                options.quality_threshold,
                options.reads_threshold,
                options.frequency_threshold)

    log.info('Done.')


if __name__ == '__main__':
    main()
