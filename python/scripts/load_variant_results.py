import csv
import sqlalchemy

# logging configuration
import log as logger

from dnascissors.config import cfg
from dnascissors.model import Base
from dnascissors.model import SequencingLibraryContent
from dnascissors.model import VariantResult


def loadVariantResults(log, session, file_name):
    log.info("Loading variant results from {:s}".format(file_name))
    with open(file_name, 'r') as f:
        reader = csv.DictReader(f, delimiter=',')
        for row in reader:
            seq_lib_content = session.query(SequencingLibraryContent) \
                                     .filter(SequencingLibraryContent.sequencing_sample_name == row['sequencing_sample_name']) \
                                     .filter(SequencingLibraryContent.sequencing_barcode == row['sequencing_barcode']) \
                                     .first()
            if not seq_lib_content:
                raise Exception("There is no sequencing library content for {:s} sample with {:s} barcode".format(row['sequencing_sample_name'], row['sequencing_barcode']))
            result = VariantResult(sequencing_library_content=seq_lib_content)
            result.variant_type = row['variant_type']
            result.gene_id = row['gene_id']
            result.cDNA_effect = row['cDNA_effect']
            result.protein_effect = row['protein_effect']
            result.codons = row['codons']
            result.chromosome = row['chromosome']
            result.position = int(row['position'])
            result.sequence_ref = row['sequence_ref']
            result.sequence_alt = row['sequence_alt']
            result.allele_fraction = float(row['allele_fraction'])
            result.depth = int(row['depth'])
            result.quality = int(row['quality'])
            result.amplicon = row['amplicon']
            result.exon = row['exon']
            result.offset_from_primer_end = int(row['offset_from_primer_end'])
            result.indel_length = int(row['indel_length'])
            result.sequence_forward_context = row['sequence_forward_context']
            result.sequence_reverse_context = row['sequence_reverse_context']
            result.sequence_alleles = row['sequence_alleles']
            result.ge_score = float(row['ge_score'])
            session.add(result)
            log.info("Variant result added for {:s} sample with {:s} barcode".format(row['sequencing_sample_name'], row['sequencing_barcode']))

    session.commit()


def main():
    import argparse
    import os

    parser = argparse.ArgumentParser()
    parser.add_argument("--file", dest="file", action="store", help="The variant results CSV file e.g. 'data/20170127_GEP00001/GEP00001_NGS_IndelsResults.csv'.", required=True)
    parser.add_argument("--clean", dest="clean_db", action="store_true", default=False, help="Clean database before loading?")
    options = parser.parse_args()

    log = logger.get_custom_logger(os.path.join(os.path.dirname(__file__), 'load_variant_results.log'))

    engine = sqlalchemy.create_engine(cfg['DATABASE_URI'])
    Base.metadata.bind = engine
    DBSession = sqlalchemy.orm.sessionmaker(bind=engine)
    session = DBSession()

    try:
        if options.clean_db:
            delete_count = session.query(VariantResult).delete()
            log.info('Deleted {:d} variant results'.format(delete_count))

        loadVariantResults(log, session, options.file)
    except Exception as e:
        log.exception(e)

    session.close()


if __name__ == '__main__':
    main()
