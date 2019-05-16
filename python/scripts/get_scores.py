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

    def get_scores(self, data_file):
        results = self.dbsession.query(Well)\
                      .join(Well.experiment_layout)\
                      .join(ExperimentLayout.project)\
                      .filter(Project.geid == self.project_geid)\
                      .all()
        # data collection
        data = []
        for i, well in enumerate(results, start=1):
            zygosity = 'empty-well'
            guide_name = 'no-guide'
            well_content_type = 'empty-well'
            score = 0
            well_abundance_ratio = 0
            if well.abundances:
                well_abundance_ratio = well.abundances[0].ratio_800_700
            if well.well_content:
                well_content_type = 'sample-well'
                zygosity = 'wt'
                if well.well_content.guides:
                    guide_name = well.well_content.guides[0].name
            if well.sequencing_library_contents:
                # TODO need a loop here to make sure dna_source is fixed cells
                # cannot trust order of items in list returned by sqlalchemy
                # sequencing_library_contents[0].dna_source corresponds to fixed cells,
                # and sequencing_library_contents[1].dna_source to gDNA
                if well.sequencing_library_contents[0].mutation_summaries:
                    score = well.sequencing_library_contents[0].mutation_summaries[0].score
                    zygosity = well.sequencing_library_contents[0].mutation_summaries[0].zygosity
            data.append([
                        well.experiment_layout.project.geid,
                        well.experiment_layout.geid,
                        well.row,
                        well.column,
                        guide_name,
                        well_abundance_ratio,
                        zygosity,
                        score,
                        well_content_type
                        ])
        if data:
            # convert to pandas dataframe and add column names
            df = pandas.DataFrame(data)
            column_names = ['project', 'plate', 'row', 'column', 'guide', 'protein_abundance', 'zygosity', 'score', 'well_content_type']
            df = df.rename(columns=dict(enumerate(column_names)))
            df.to_csv(data_file, index=False)


def main():
    logging.basicConfig(level=logging.DEBUG,
                        format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                        datefmt='%m-%d %H:%M')
    engine = sqlalchemy.create_engine(cfg['DATABASE_URI'])
    Base.metadata.bind = engine
    DBSession = sqlalchemy.orm.sessionmaker(bind=engine)
    session = DBSession()
    try:
        data = DataExtractor(session, 'GEP00001')
        data.get_scores("scores.csv")

    except Exception as e:
        logging.exception(e)
    finally:
        session.close()

if __name__ == '__main__':
    main()
