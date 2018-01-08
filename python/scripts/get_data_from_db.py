import sqlalchemy
import json

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
                 .outerjoin(SequencingLibraryContent.variant_results)\
                 .join(SequencingLibraryContent.well)\
                 .join(Well.well_content)\
                 .join(Well.experiment_layout)\
                 .join(ExperimentLayout.project)\
                 .filter(Project.geid == 'GEP00001')\
                 .all()

# Growth data: Unique ID | hours | confluence_percentage
# Protein data: Unique ID | intensity_channel_700 | intensity_channel_800
# Variant results: Unique ID | variant_caller | variant_type | chromosome | position | ref | alt | allele_fraction | depth | amplicon | indel_length
# Mutation summaries: Unique ID | zygosity | consequence

data = []
for sequencing_library_content in results:
    well = sequencing_library_content.well
    row = {'well_id': well.id,
           'well_location': '{}_{}{}'.format(well.experiment_layout.geid, well.row, well.column),
           'sample_name': sequencing_library_content.sequencing_sample_name,
           'abundances': [{'700': a.intensity_channel_700, '800': a.intensity_channel_800} for a in well.abundances],
           'growths': [(g.hours, g.confluence_percentage) for g in well.growths],
           'variant_results': [(v.variant_caller, v.variant_type, v.chromosome, v.position, v.ref, v.alt, v.allele_fraction, v.depth, v.amplicon, v.indel_length) for v in sequencing_library_content.variant_results],
           'mutation_summaries': [(m.zygosity, m.consequence, m.variant_caller_presence, m.score) for m in sequencing_library_content.mutation_summaries]
          }
    data.append(row)

with open('data_from_db.json', 'w') as f:
    for d in data:
        f.write(json.dumps(d, indent=2))
