# -*- coding: utf-8 -*-
"""
Created on Tue Jul 11 12:32:07 2017

@author: alvare01

plot description: 
Query database to get sample names, plate, well locations, guides, scores, protein value, slope value and zygosity per clone
Create a scatter plot with the shape of a 96-well plate, mapping the scores to a colour gradient

"""

import sqlalchemy
import pandas
import plotly.graph_objs as go
import plotly.offline as pyoff

from dnascissors.config import cfg
from dnascissors.model import Base
from dnascissors.model import SequencingLibraryContent
from dnascissors.model import Well
from dnascissors.model import ExperimentLayout
from dnascissors.model import Project

#-------- sqlalchemy query
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
                 .filter(SequencingLibraryContent.dna_source != 'gDNA')\
                 .all()

data = []
zygosities = []
for i in results:
    mut_zygosities = 'wt'
    guide_name = 'none'
    well = i.well
    row = [i.sequencing_sample_name,
           i.well.experiment_layout.geid,
           well.row,
           well.column,
           ':'.join([g.name for g in well.well_content.guides]),
            well.abundances[0].intensity_channel_800/well.abundances[0].intensity_channel_700
           #[a.intensity_channel_800/a.intensity_channel_700 for a in well.abundances],
           #[(g.hours, g.confluence_percentage) for g in well.growths],
           ]
    if i.mutation_summaries:
        mut_zygosities = i.mutation_summaries[0].zygosity
    rowzygo = row + [mut_zygosities]
    data.append(rowzygo)

data = pandas.DataFrame(data)
data.columns['sample', 'plate', 'row', 'column', 'guide', 'protein', 'zygosity']

