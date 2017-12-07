"""
Python3 file of the genome-editing project
Created by Anne Pajon @pajanne on 08/03/2017
"""

from dnascissors.model import Base
from dnascissors.config import cfg
from sqlalchemy import create_engine

engine = create_engine(cfg['DATABASE_URI'])
Base.metadata.drop_all(engine)
Base.metadata.create_all(engine)
