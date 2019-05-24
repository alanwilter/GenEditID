import pandas
import sqlalchemy
import logging

from dnascissors.config import cfg
from dnascissors.model import Base
from dnascissors.model import ExperimentLayout
from dnascissors.model import Project
from dnascissors.model import Well


class DataExtractor:

    def __init__(self, dbsession, project_geid):
        self.include_js = False
        self.dbsession = dbsession
        self.project_geid = project_geid

    def get_data(self, data_file):
        results = self.dbsession.query(Well)\
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
            column_names = ['plate_id', 'well', 'cell_line_name', 'clone_name', 'guide_name', 'is_control', 'content_type', 'dna_source', 'sequencing_barcode', 'intensity_channel_700', 'intensity_channel_800', 'protein_abundance']
            df = df.rename(columns=dict(enumerate(column_names)))
            df.to_csv(data_file, index=False)

        if barcodes:
            # convert to pandas dataframe and add column names
            df = pandas.DataFrame(barcodes)
            column_names = ['plate_id', 'well', 'sequencing_barcode']
            df = df.rename(columns=dict(enumerate(column_names)))
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
        data = DataExtractor(session, 'GEP00005')
        data.get_data("expression_data.csv")

    except Exception as e:
        logging.exception(e)
    finally:
        session.close()

if __name__ == '__main__':
    main()
