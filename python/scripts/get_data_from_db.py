import sqlalchemy

from dnascissors.config import cfg
from dnascissors.model import Base
from dnascissors.model import SequencingLibraryContent
from dnascissors.model import Well
from dnascissors.model import ExperimentLayout
from dnascissors.model import Project

# sqlalchemy database connection
engine = sqlalchemy.create_engine(cfg['DATABASE_URI'])
Base.metadata.bind = engine
DBSession = sqlalchemy.orm.sessionmaker(bind=engine)
session = DBSession()

# sqlalchemy query
results = session.query(SequencingLibraryContent)\
                 .outerjoin(SequencingLibraryContent.mutation_summaries)\
                 .join(SequencingLibraryContent.well)\
                 .join(Well.well_content)\
                 .join(Well.experiment_layout)\
                 .join(ExperimentLayout.project)\
                 .filter(Project.geid == 'GEP00001')\
                 .all()

data = []
for sequencing_library_content in results:
    well = sequencing_library_content.well
    row = [sequencing_library_content.sequencing_sample_name, \
           sequencing_library_content.well.experiment_layout.geid, \
           '{}{}'.format(well.row, well.column), \
           ':'.join([g.name for g in well.well_content.guides]), \
           [a.intensity_channel_700 for a in well.abundances], \
           [a.intensity_channel_800 for a in well.abundances], \
           [(g.hours, g.confluence_percentage) for g in well.growths], \
           [m.zygosity for m in sequencing_library_content.mutation_summaries]
           ]
    data.append(row)

for d in data:
    print(d)
