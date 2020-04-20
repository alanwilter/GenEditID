import os
import glob
import pandas
import sqlalchemy
import logging
import sys

from geneditid.config import cfg
from geneditid.model import Base
from geneditid.model import ExperimentLayout
from geneditid.model import Project
from geneditid.model import Well


class DataExtractor:

    def __init__(self, dbsession, project_geid):
        self.dbsession = dbsession
        self.project_geid = project_geid

    def get_data(self, data_file):
        results = self.dbsession.query(Well)\
                      .join(Well.well_content)\
                      .join(Well.experiment_layout)\
                      .join(ExperimentLayout.project)\
                      .filter(Project.geid == self.project_geid)\
                      .all()
        # data collection
        data = []
        barcodes = []
        for i, well in enumerate(results, start=1):
            well_content_type = 'empty-well'
            is_control = False
            guide_name = 'no-guide'
            clone_name = 'no-clone'
            cell_line_name = 'no-cell-line'
            dna_source = 'no-dna-source'
            sequencing_barcode = 'no-barcode'
            intensity_channel_700 = 0
            intensity_channel_800 = 0
            well_abundance_ratio = 0
            if well.abundances:
                intensity_channel_700 = well.abundances[0].intensity_channel_700
                intensity_channel_800 = well.abundances[0].intensity_channel_800
                well_abundance_ratio = well.abundances[0].ratio_800_700
            if well.well_content:
                well_content_type = well.well_content.content_type
                is_control = well.well_content.is_control
                if well.well_content.guides:
                    guide_name = well.well_content.guides[0].name
                if well.well_content.clone:
                    clone_name = well.well_content.clone.name
                    if well.well_content.clone.cell_line:
                        cell_line_name = well.well_content.clone.cell_line.name
            if well.sequencing_library_contents:
                dna_source = well.sequencing_library_contents[0].dna_source
                sequencing_barcode = well.sequencing_library_contents[0].sequencing_barcode
            data.append([
                        well.experiment_layout.geid,
                        '{}{}'.format(well.row, well.column),
                        cell_line_name,
                        clone_name,
                        guide_name,
                        is_control,
                        well_content_type,
                        dna_source,
                        sequencing_barcode,
                        intensity_channel_700,
                        intensity_channel_800,
                        well_abundance_ratio
                        ])
            barcodes.append([well.experiment_layout.geid,
                             '{}{}'.format(well.row, well.column),
                             sequencing_barcode])
        if data:
            # convert to pandas dataframe and add column names
            df = pandas.DataFrame(data)
            column_names = ['plate_id', 'well', 'cell_line_name', 'clone_name', 'guide_name', 'is_control', 'content_type', 'dna_source', 'sample_id', 'intensity_channel_700', 'intensity_channel_800', 'protein_abundance']
            df = df.rename(columns=dict(enumerate(column_names)))
            df = df[df['protein_abundance'] != 0]
            if not df.empty:
                df.to_csv(data_file, index=False)

        if barcodes:
            # convert to pandas dataframe and add column names
            df = pandas.DataFrame(barcodes)
            column_names = ['plate_id', 'well', 'sequencing_barcode']
            df = df.rename(columns=dict(enumerate(column_names)))
            if self.project_geid == 'GEP00005':
                df = df[df['plate_id'] == 'GEP00005_02']
            df.to_csv('barcodes.csv', index=False)


def main():
    logging.basicConfig(level=logging.DEBUG,
                        format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                        datefmt='%m-%d %H:%M')
    engine = sqlalchemy.create_engine(cfg['DATABASE_URI'])
    Base.metadata.bind = engine
    DBSession = sqlalchemy.orm.sessionmaker(bind=engine)
    session = DBSession()
    try:
        geid = sys.argv[1]  # add your GEPID on the command line
        data = DataExtractor(session, geid)
        data.get_data("expression_data.csv")

        varianid_folder = 'editid_variantid'
        df_barcodes = pandas.read_csv('barcodes.csv')

        df_variants = pandas.read_csv(os.path.join(varianid_folder, 'variantid.csv'))
        df_impacts = pandas.read_csv(os.path.join(varianid_folder, 'impacts.csv'))

        df_variants = df_variants.merge(df_barcodes, left_on='sample_id', right_on='sequencing_barcode', how='left')
        df_impacts = df_impacts.merge(df_barcodes, left_on='sample_id', right_on='sequencing_barcode', how='left')

        df_variants = df_variants[['plate_id', 'well', 'sample_id', 'amplicon_id', 'total_reads', 'amplicon_reads', 'amplicon_filtered_reads', 'amplicon_low_quality_reads', 'amplicon_primer_dimer_reads', 'amplicon_low_abundance_reads', 'variant_reads', 'variant_frequency', 'sequence', 'variant_id', 'variant_type', 'variant_consequence', 'variant_score']]
        df_variants.to_csv('editid_variantid/variantid_with_plate_location.csv', index=False)

        df_impacts = df_impacts[['plate_id', 'well', 'sample_id', 'amplicon_id', 'impact', 'impact_frequency']]
        df_impacts.to_csv('editid_variantid/impacts_with_plate_location.csv', index=False)

        for file in glob.glob(os.path.join(varianid_folder, 'koscores_*.csv')):
            if 'with_plate_location' not in file:
                df_koscores = pandas.read_csv(file)
                df_koscores = df_koscores.merge(df_barcodes, left_on='sample_id', right_on='sequencing_barcode', how='left')
                df_koscores = df_koscores[['plate_id', 'well', 'sample_id', 'HighImpact', 'MediumImpact', 'LowImpact', 'WildType', 'LowFrequency', 'koscore']]
                output, ext = os.path.splitext(os.path.basename(file))
                df_koscores.to_csv(os.path.join(varianid_folder, '{}_with_plate_location{}'.format(output, ext)), index=False)

    except Exception as e:
        logging.exception(e)
    finally:
        session.close()

if __name__ == '__main__':
    main()
