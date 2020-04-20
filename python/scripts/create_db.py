"""
Python3 file of the genome-editing project
Created by Anne Pajon @pajanne on 08/03/2017
"""
import os
from geneditid.model import Base
from geneditid.config import cfg
from sqlalchemy import create_engine

engine = create_engine(cfg['DATABASE_URI'])
Base.metadata.drop_all(engine)
Base.metadata.create_all(engine)

# create projects folder
if not os.path.exists(cfg['PROJECTS_FOLDER']):
    os.makedirs(cfg['PROJECTS_FOLDER'])
