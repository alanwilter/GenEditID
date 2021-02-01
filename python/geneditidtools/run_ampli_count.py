import os
import gzip
import argparse
import collections

import pandas as pd
from Bio import SeqIO

import geneditid.log as logger
from geneditid.config import cfg


def count_reads(log, outputfile, fastq_dir, fastq_extension, amplicons, quality_threshold, abundance_threshold, reverse_flag, targets=pd.DataFrame()):
    filename, ext = os.path.splitext(outputfile)
    outputfile_tsearch = '{}_tsearch{}'.format(filename, ext)
    with open(outputfile, 'w') as out, open(outputfile_tsearch, 'w') as out_tsearch:
        out.write("sample_id,amplicon_id,total_reads,amplicon_reads,amplicon_filtered_reads,amplicon_low_quality_reads,amplicon_primer_dimer_reads,amplicon_low_abundance_reads,variant_reads,variant_frequency,sequence\n")
        out_tsearch.write("sample_id,amplicon_id,total_reads,amplicon_reads,amplicon_filtered_reads,amplicon_low_quality_reads,amplicon_primer_dimer_reads,amplicon_low_abundance_reads,variant_reads,variant_frequency,sequence,tsearch_id\n")
        for filename in sorted(os.listdir(fastq_dir)):
            if filename.endswith(fastq_extension):
                splited_filename = filename.split('.')
                slxid, sample_id = splited_filename[0], splited_filename[1]
                total_reads = 0
                amplicon_reads = collections.OrderedDict()
                amplicon_filtered_reads = collections.OrderedDict()
                amplicon_low_quality_reads = collections.OrderedDict()
                amplicon_primer_dimer_reads = collections.OrderedDict()
                amplicon_low_abundance_reads = collections.OrderedDict()
                variant_reads = collections.OrderedDict()
                for i, amplicon in amplicons.iterrows():
                    variant_reads[amplicon['id']] = {}
                    amplicon_filtered_reads[amplicon['id']] = 0
                    amplicon_low_quality_reads[amplicon['id']] = 0
                    amplicon_primer_dimer_reads[amplicon['id']] = 0
                    amplicon_low_abundance_reads[amplicon['id']] = 0
                # classifying and filtering reads
                with gzip.open(os.path.join(fastq_dir, filename), "rt") as handle:
                    log.info('Counting reads in file {}'.format(filename))
                    for r in SeqIO.parse(handle, "fastq"):
                        if reverse_flag:
                            rec = r.reverse_complement(id="REVCOMP")
                        else:
                            rec = r
                        total_reads += 1
                        for i, amplicon in amplicons.iterrows():
                            fprimer_pos = str(rec.seq).find(amplicon['fprimer'])
                            rprimer_pos = str(rec.seq).find(amplicon['rprimer'])
                            # classify reads per amplicon based on forward and reverse primer sequences
                            if not fprimer_pos == -1 and not rprimer_pos == -1:
                                amplicon_reads[amplicon['id']] = amplicon_reads.get(amplicon['id'], 0) + 1
                                rprimer_pos = rprimer_pos + len(amplicon['rprimer'])
                                variant_seq = str(rec.seq)[fprimer_pos:rprimer_pos]
                                # count reads associated to primer-dimer
                                if len(variant_seq) > len(amplicon['fprimer']) + len(amplicon['rprimer']) + 10:
                                    # count reads with low-quality over a window of 5 bases above the quality_threshold
                                    qualities = pd.DataFrame(data={'quality': rec.letter_annotations["phred_quality"][fprimer_pos:rprimer_pos]})
                                    rolling_quality_means = qualities.rolling(window=5, min_periods=1, on='quality').mean()
                                    if min(rolling_quality_means['quality']) >= int(quality_threshold):
                                        variant_reads[amplicon['id']][variant_seq] = variant_reads[amplicon['id']].get(variant_seq, 0) + 1
                                        amplicon_filtered_reads[amplicon['id']] = amplicon_filtered_reads.get(amplicon['id'], 0) + 1
                                    else:
                                        amplicon_low_quality_reads[amplicon['id']] = amplicon_low_quality_reads.get(amplicon['id'], 0) + 1
                                else:
                                    amplicon_primer_dimer_reads[amplicon['id']] = amplicon_primer_dimer_reads.get(amplicon['id'], 0) + 1

                    log.info("{},{},{}".format(sample_id, total_reads, ','.join(str(n) for n in amplicon_filtered_reads.values())))

                for amplicon_id in variant_reads.keys():
                    for seq in variant_reads[amplicon_id].keys():
                        # count variant reads with low abundance
                        if variant_reads[amplicon_id][seq] <= abundance_threshold:
                            amplicon_low_abundance_reads[amplicon['id']] = amplicon_low_abundance_reads.get(amplicon['id'], 0) + variant_reads[amplicon_id][seq]

                # write outputs
                for amplicon_id in variant_reads.keys():
                    # write 'amplicount.csv' if above a certain abundance threshold
                    for seq in variant_reads[amplicon_id].keys():
                        if variant_reads[amplicon_id][seq] > abundance_threshold:
                            out.write("{},{},{},{},{},{},{},{},{},{:.2f},{}\n".format(sample_id,
                                                                                      amplicon_id,
                                                                                      total_reads,
                                                                                      amplicon_reads[amplicon_id],
                                                                                      amplicon_filtered_reads[amplicon_id],
                                                                                      amplicon_low_quality_reads[amplicon_id],
                                                                                      amplicon_primer_dimer_reads[amplicon_id],
                                                                                      amplicon_low_abundance_reads[amplicon_id],
                                                                                      variant_reads[amplicon_id][seq],
                                                                                      (variant_reads[amplicon_id][seq] * 100) / amplicon_filtered_reads[amplicon_id],
                                                                                      seq))
                    # write 'amplicount_tsearch.csv' if sequence in list of desired sequences
                    if not targets.empty:
                        for i, target in targets.iterrows():
                            if target['amplicon_id'] == amplicon_id:
                                for seq in variant_reads[amplicon_id].keys():
                                    if seq.count(target['sequence']) > 0:
                                        out_tsearch.write("{},{},{},{},{},{},{},{},{},{:.2f},{},{}\n".format(sample_id,
                                                                                                             amplicon_id,
                                                                                                             total_reads,
                                                                                                             amplicon_reads[amplicon_id],
                                                                                                             amplicon_filtered_reads[amplicon_id],
                                                                                                             amplicon_low_quality_reads[amplicon_id],
                                                                                                             amplicon_primer_dimer_reads[amplicon_id],
                                                                                                             amplicon_low_abundance_reads[amplicon_id],
                                                                                                             variant_reads[amplicon_id][seq],
                                                                                                             (variant_reads[amplicon_id][seq] * 100) / amplicon_filtered_reads[amplicon_id],
                                                                                                             seq,
                                                                                                             target['id']))


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("--config", dest="config", action="store", help="The 4 columns input config file: 'id,fprimer,rprimer,amplicon'", default='amplicount_config.csv', required=False)
    parser.add_argument("--fastqdir", dest="fastq_dir", action="store", help="Fastq file directory", default='fastq', required=False)
    parser.add_argument("--fastqext", dest="fastq_extension", action="store", help="Fastq file extension", default='.fqjoin.gz', required=False)
    parser.add_argument("--quality", dest="quality_threshold", action="store", help="Quality threshold for average phred quality across a window over the amplicon sequence", default=10, required=False)
    parser.add_argument("--abundance", dest="abundance_threshold", action="store", help="Abundance threshold for min number of reads to report per variant", default=60, required=False)
    parser.add_argument("--reverse", dest="reverse_complement", action="store_true", help="Reverse complement all reads", required=False)
    parser.add_argument("--output", dest="output", action="store", help="The output file", default='amplicount.csv', required=False)
    parser.add_argument("--withseq", dest="sequences", action="store", help="The 3 columns input sequence file: 'amplicon_id,id,sequence'", default='amplicount_config_tsearch.csv', required=False)
    options = parser.parse_args()

    log = logger.get_custom_logger('amplicount.log')

    log.info('>>> Getting list of amplicons...')
    amplicons = pd.read_csv(options.config)
    if amplicons.empty:
        log.warning('No amplicon found in file {}'.format(options.config))
    else:
        for i, amplicon in amplicons.iterrows():
            log.info('Amplicon: {}'.format(amplicon['id']))
            log.info('Forward primer: {} {}'.format(len(amplicon['fprimer']), amplicon['fprimer']))
            log.info('Reverse primer: {} {}'.format(len(amplicon['rprimer']), amplicon['rprimer']))
            log.info('Ref amplicon seq: {} {}'.format(len(amplicon['amplicon']), amplicon['amplicon']))
            log.info('----------')

    targets = pd.DataFrame()
    if os.path.exists(options.sequences):
        log.info('>>> Getting list of desire edited sequences...')
        targets = pd.read_csv(options.sequences)
        if targets.empty:
            log.warning('No target found in file {}'.format(options.sequences))
        else:
            for i, target in targets.iterrows():
                log.info('Desire edited sequence: {} {}'.format(target['amplicon_id'], target['id'], target['sequence']))

    log.info('>>> Counting reads...')
    log.info('FastQ file directory: {}'.format(options.fastq_dir))
    log.info('FastQ extension: {}'.format(options.fastq_extension))
    log.info('Read quality Threshold: {}'.format(options.quality_threshold))
    log.info('Abundance threshold: {}'.format(options.abundance_threshold))
    log.info('Count reads with a minimum base quality score above the read quality threshold.')
    log.info('Report only variants with variant read counts which is above the aboundance threshold.')
    log.info('Report also variants with sequence that matches the desire edited sequence if provided.')

    count_reads(log,
                options.output,
                options.fastq_dir,
                options.fastq_extension,
                amplicons,
                float(options.quality_threshold),
                int(options.abundance_threshold),
                options.reverse_complement,
                targets)

    log.info('Done.')

if __name__ == '__main__':
    main()
