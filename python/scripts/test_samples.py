import sqlalchemy

from dnascissors.config import cfg
from dnascissors.model import Base
from dnascissors.model import Project
from dnascissors.model import ExperimentLayout
from dnascissors.model import Well
from dnascissors.model import SequencingLibraryContent

engine = sqlalchemy.create_engine(cfg['DATABASE_URI'])
Base.metadata.bind = engine
DBSession = sqlalchemy.orm.sessionmaker(bind=engine)
dbsession = DBSession()
geproject = dbsession.query(Project).filter(Project.geid == "GEP00011").one()
samples = dbsession\
              .query(SequencingLibraryContent)\
              .join(SequencingLibraryContent.well)\
              .join(Well.experiment_layout)\
              .join(ExperimentLayout.project)\
              .filter(Project.geid == 'GEP00011')\
              .all()

for s in samples:
    print(s.dna_source)
