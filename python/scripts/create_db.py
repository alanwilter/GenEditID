import os
from sqlalchemy import create_engine
from geneditid.model import Base
from geneditid.config import cfg

engine = create_engine(cfg['DATABASE_URI'])
Base.metadata.drop_all(engine)
Base.metadata.create_all(engine)

# create projects folder
if not os.path.exists(cfg['PROJECTS_FOLDER']):
    os.makedirs(cfg['PROJECTS_FOLDER'])
