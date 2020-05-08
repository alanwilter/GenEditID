import sqlalchemy
from geneditid.model import Base
from geneditid.config import cfg

engine = sqlalchemy.create_engine(cfg['DATABASE_URI'])
Base.metadata.bind = engine
DBSession = sqlalchemy.orm.sessionmaker(bind=engine)
dbsession = DBSession()
