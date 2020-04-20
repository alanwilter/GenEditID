import os
import gzip
import argparse
import collections
import pandas as pd
import scripts.log as logger
from Bio import SeqIO
from geneditid.config import cfg


def count_reads(log, outputfile, fastq_dir, fastq_extension, amplicons, quality_threshold, abundance_threshold):
    filename, ext = os.path.splitext(outputfile)
    with open(outputfile, 'w') as out:
        out.write("sample_id,amplicon_id,total_reads,amplicon_reads,amplicon_filtered_reads,amplicon_low_quality_reads,amplicon_primer_dimer_reads,amplicon_low_abundance_reads,variant_reads,variant_frequency,sequence\n")
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
                    for rec in SeqIO.parse(handle, "fastq"):
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
                                    #rolling_quality_means.dropna(inplace=True)
                                    #log.debug(rolling_quality_means)
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

                # write output
                for amplicon_id in variant_reads.keys():
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


def count_reads_for_target_sequences(log, outputfile, fastq_dir, fastq_extension, targets, quality_threshold):
    filename, ext = os.path.splitext(outputfile)
    with open(outputfile, 'w') as out:
        out.write("sample_id,target_id,total_reads,target_reads,target_filtered_reads,target_low_quality_reads,sequence\n")
        for filename in sorted(os.listdir(fastq_dir)):
            if filename.endswith(fastq_extension):
                splited_filename = filename.split('.')
                slxid, sample_id = splited_filename[0], splited_filename[1]
                total_reads = 0
                target_reads = collections.OrderedDict()
                target_filtered_reads = collections.OrderedDict()
                target_low_quality_reads = collections.OrderedDict()
                for i, target in targets.iterrows():
                    target_reads[target['id']] = 0
                    target_filtered_reads[target['id']] = 0
                    target_low_quality_reads[target['id']] = 0
                # classifying and filtering reads
                with gzip.open(os.path.join(fastq_dir, filename), "rt") as handle:
                    log.info('Counting reads in file {}'.format(filename))
                    for rec in SeqIO.parse(handle, "fastq"):
                        total_reads += 1
                        for i, target in targets.iterrows():
                            target_pos = str(rec.seq).find(target['sequence'])
                            if not target_pos == -1:
                                target_reads[target['id']] = target_reads.get(target['id'], 0) + 1
                                # count reads with low-quality over a window of 5 bases above the quality_threshold
                                qualities = pd.DataFrame(data={'quality': rec.letter_annotations["phred_quality"][target_pos:target_pos+len(target['sequence'])]})
                                rolling_quality_means = qualities.rolling(window=5, min_periods=1, on='quality').mean()
                                #rolling_quality_means.dropna(inplace=True)
                                #log.debug(rolling_quality_means)
                                if min(rolling_quality_means['quality']) >= int(quality_threshold):
                                    target_filtered_reads[target['id']] = target_filtered_reads.get(target['id'], 0) + 1
                                else:
                                    target_low_quality_reads[target['id']] = target_low_quality_reads.get(target['id'], 0) + 1

                    log.info("{},{},{}".format(sample_id, total_reads, ','.join(str(n) for n in target_filtered_reads.values())))

                # write output
                for target_id in target_reads.keys():
                    for i, target in targets.iterrows():
                        if target_reads[target_id] > 0 and target_id == target['id']:
                            out.write("{},{},{},{},{},{},{}\n".format(sample_id,
                                                                      target_id,
                                                                      total_reads,
                                                                      target_reads[target_id],
                                                                      target_filtered_reads[target_id],
                                                                      target_low_quality_reads[target_id],
                                                                      target['sequence']))


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("--config", dest="config", action="store", help="The 4 columns input config file: 'id,fprimer,rprimer,amplicon'", default='amplicount_config.csv', required=False)
    parser.add_argument("--fastqdir", dest="fastq_dir", action="store", help="Fastq file directory", required=True)
    parser.add_argument("--fastqext", dest="fastq_extension", action="store", help="Fastq file extension", default='.fqjoin.gz', required=False)
    parser.add_argument("--quality", dest="quality_threshold", action="store", help="Quality threshold for average phred quality across a window over the amplicon sequence", default=10, required=False)
    parser.add_argument("--abundance", dest="abundance_threshold", action="store", help="Abundance threshold for min number of reads to report per variant", default=60, required=False)
    parser.add_argument("--output", dest="output", action="store", help="The output file", default='amplicount.csv', required=False)
    parser.add_argument("--with_seq", dest="sequences", action="store", help="The 2 columns input sequence file: 'id,sequence'", required=False)
    options = parser.parse_args()

    log = logger.get_custom_logger(os.path.join(cfg['PROJECTS_FOLDER'], 'amplicount.log'))

    if not options.sequences:
        log.info('>>> Getting list of amplicons...')
        amplicons = pd.read_csv(options.config)
        for i, amplicon in amplicons.iterrows():
            log.info('Amplicon: {}'.format(amplicon['id']))
            log.info('Forward primer: {} {}'.format(len(amplicon['fprimer']), amplicon['fprimer']))
            log.info('Reverse primer: {} {}'.format(len(amplicon['rprimer']), amplicon['rprimer']))
            log.info('Ref amplicon seq: {} {}'.format(len(amplicon['amplicon']), amplicon['amplicon']))
            log.info('----------')

        log.info('>>> Counting reads per variant per amplicon...')
        log.info('FastQ file directory: {}'.format(options.fastq_dir))
        log.info('FastQ extension: {}'.format(options.fastq_extension))
        log.info('Read quality Threshold: {}'.format(options.quality_threshold))
        log.info('Read abundance Threshold: {}'.format(options.abundance_threshold))
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
                    options.abundance_threshold)

        log.info('Done.')

    else:
        # count reads for target sequences provided
        log.info('>>> Getting list of target sequences...')
        targets = pd.read_csv(options.sequences)
        for i, target in targets.iterrows():
            log.info('Target sequence: {} {}'.format(target['id'], target['sequence']))

        log.info('>>> Counting reads per target sequences...')
        log.info('FastQ file directory: {}'.format(options.fastq_dir))
        log.info('FastQ extension: {}'.format(options.fastq_extension))
        log.info('Read quality Threshold: {}'.format(options.quality_threshold))
        # count reads associated with each target sequence
        count_reads_for_target_sequences(log,
                                         options.output,
                                         options.fastq_dir,
                                         options.fastq_extension,
                                         targets,
                                         options.quality_threshold)

        log.info('Done.')

if __name__ == '__main__':
    main()
