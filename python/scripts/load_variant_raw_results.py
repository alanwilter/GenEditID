import csv
import sqlalchemy
import pandas

# logging configuration
import log as logger

from dnascissors.config import cfg
from dnascissors.model import Base
from dnascissors.model import SequencingLibraryContent
from dnascissors.model import VariantRawResult
from dnascissors.model import VariantResult


def to_int(text):
    if pandas.isnull(text):
        return None
    else:
        try:
            return int(text)
        except ValueError:
            return None


def load_variant_raw_results_per_sheet(log, sheet, session, variant_caller, variant_type):
    header = ['sample',
              'barcode',
              'variant_type_consequence',
              'symbol_gene_id',
              'cdna_effect',
              'protein_effect',
              'codons',
              'chromosome',
              'position',
              'ref',
              'alt',
              'igv_links',
              'specific',
              'confidence',
              'allele_fraction',
              'depth',
              'filter',
              'quality',
              'amplicon',
              'pubmed',
              'gene',
              'exon',
              'intron',
              'existing_variation',
              'sift',
              'polyphen',
              'clinical_significance',
              'gmaf',
              'offset_from_primer_end',
              'indel_length',
              'forward_context',
              'alleles',
              'reverse_context']
    sheet.columns = header
    for i, row in enumerate(sheet.itertuples(), 1):
        log.debug(row)
        seq_lib_content = session.query(SequencingLibraryContent) \
                                 .filter(SequencingLibraryContent.sequencing_sample_name == row.sample) \
                                 .filter(SequencingLibraryContent.sequencing_barcode == row.barcode) \
                                 .first()
        if not seq_lib_content:
            raise Exception("There is no sequencing library content for {:s} sample with {:s} barcode".format(row.sample, row.barcode))
        result = VariantRawResult(sequencing_library_content=seq_lib_content)
        result.variant_caller = variant_caller
        result.variant_type = variant_type
        result.consequence = row.variant_type_consequence
        result.gene_id = row.symbol_gene_id
        result.cdna_effect = row.cdna_effect
        result.protein_effect = row.protein_effect
        result.codons = row.codons
        result.chromosome = row.chromosome
        result.position = to_int(row.position)
        result.ref = row.ref
        result.alt = row.alt
        result.allele_fraction = float(row.allele_fraction)
        result.depth = to_int(row.depth)
        result.quality = to_int(row.quality)
        result.amplicon = row.amplicon
        result.gene = row.gene
        result.exon = row.exon
        result.intron = row.intron
        result.existing_variation = row.existing_variation
        result.sift = row.sift
        result.polyphen = row.polyphen
        result.clinical_significance = row.clinical_significance
        result.gmaf = row.gmaf
        result.offset_from_primer_end = to_int(row.offset_from_primer_end)
        result.indel_length = to_int(row.indel_length)
        result.forward_context = row.forward_context
        result.alleles = row.alleles
        result.reverse_context = row.reverse_context
        session.add(result)
        log.info("Variant result added for {:s} sample with {:s} barcode".format(row.sample, row.barcode))
    session.commit()


def load_variant_raw_results(log, session, workbook_file, variant_caller):
    xls = pandas.ExcelFile(workbook_file)
    load_variant_raw_results_per_sheet(log, xls.parse('SNVs', header=None, skiprows=1), session, variant_caller, 'SNV')
    load_variant_raw_results_per_sheet(log, xls.parse('Indels', header=None, skiprows=1), session, variant_caller, 'INDEL')


def main():
    import argparse
    import os

    parser = argparse.ArgumentParser()
    parser.add_argument("--file", dest="file", action="store", help="The variant raw results excel file e.g. 'data/20170127_GEP00001/SLX-13775.vardict.variants.xlsx'.", required=True)
    parser.add_argument("--caller", dest="variant_caller", action="store", help="The type of variant caller used.", choices=['VarDict', 'HaplotypeCaller', 'FreeBayes'], required=True)
    parser.add_argument("--clean", dest="clean_db", action="store_true", default=False, help="Clean database before loading?")
    options = parser.parse_args()

    log = logger.get_custom_logger(os.path.join(os.path.dirname(__file__), 'load_variant_results.log'))

    engine = sqlalchemy.create_engine(cfg['DATABASE_URI'])
    Base.metadata.bind = engine
    DBSession = sqlalchemy.orm.sessionmaker(bind=engine)
    session = DBSession()

    try:
        if options.clean_db:
            delete_count = session.query(VariantRawResult).delete()
            log.info('Deleted {:d} variant results'.format(delete_count))

        #load_variant_results(log, session, options.file, options.variant_type)
        load_variant_raw_results(log, session, options.file, options.variant_caller)
    except Exception as e:
        log.exception(e)

    session.close()


if __name__ == '__main__':
    main()
