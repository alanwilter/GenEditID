import csv
import sqlalchemy

# logging configuration
import log as logger

from dnascissors.config import cfg
from dnascissors.model import Base
from dnascissors.model import SequencingLibraryContent
from dnascissors.model import VariantResult


def load_variant_results(log, session, file_name, variant_type):
    log.info("Loading variant results from {:s}".format(file_name))
    with open(file_name, 'r') as f:
        reader = csv.DictReader(f, delimiter=',')
        for row in reader:
            seq_lib_content = session.query(SequencingLibraryContent) \
                                     .filter(SequencingLibraryContent.sequencing_sample_name == row['Sample']) \
                                     .filter(SequencingLibraryContent.sequencing_barcode == row['Barcode']) \
                                     .first()
            if not seq_lib_content:
                raise Exception("There is no sequencing library content for {:s} sample with {:s} barcode".format(row['sequencing_sample_name'], row['sequencing_barcode']))
            result = VariantResult(sequencing_library_content=seq_lib_content)
            result.variant_type = variant_type
            result.consequence = row['Variant_type_consequence']
            result.gene_id = row['Symbol_Gene_ID']
            result.gene = row['Gene']
            result.cdna_effect = row['cDNA_effect']
            result.protein_effect = row['Protein_effect']
            result.codons = row['Codons']
            result.chromosome = row['Chromosome']
            result.position = int(row['Position'])
            result.sequence_ref = row['Ref']
            result.sequence_alt = row['Alt']
            result.allele_fraction = float(row['Allele_fraction'])
            result.depth = int(row['Depth'])
            result.quality = int(row['Quality'])
            result.amplicon = row['Amplicon']
            result.exon = row['Exon']
            result.offset_from_primer_end = int(row['Offset_from_primer_end'])
            result.indel_length = None
            if variant_type == 'INDEL':
                result.indel_length = int(row['Indel_length'])
            result.sequence_forward_context = row['fivePrimeContext']
            result.sequence_reverse_context = row['threePrimeContext']
            result.sequence_alleles = row['Alleles']
            result.on_target = bool(row['OnTarget'])
            result.cut_site = int(row['cutSite'])
            result.distance = int(row['distance'])
            result.ge_score = int(row['disScore'])
            session.add(result)
            log.info("Variant result added for {:s} sample with {:s} barcode".format(row['Sample'], row['Barcode']))

    session.commit()


def main():
    import argparse
    import os

    parser = argparse.ArgumentParser()
    parser.add_argument("--file", dest="file", action="store", help="The variant results CSV file e.g. 'data/20170127_GEP00001/GEP00001_NGS_IndelsResults.csv'.", required=True)
    parser.add_argument("--type", dest="variant_type", action="store", help="The type of variant results.", choices=['INDEL', 'SNV'], required=True)
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

        load_variant_results(log, session, options.file, options.variant_type)
    except Exception as e:
        log.exception(e)

    session.close()


if __name__ == '__main__':
    main()
