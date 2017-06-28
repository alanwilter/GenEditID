# -*- coding: utf-8 -*-
"""
Created on Wed Jun 28 14:27:18 2017

@author: alvare01
"""

import sqlalchemy
import plotly
import plotly.graph_objs as go
import plotly.offline as pyoff

from dnascissors.model import MutationSummary
from dnascissors.model import Base
from dnascissors.config import cfg
from dnascissors.model import SequencingLibraryContent
from dnascissors.model import Well
from dnascissors.model import ExperimentLayout
from dnascissors.model import Project


engine = sqlalchemy.create_engine(cfg['DATABASE_URI'])
Base.metadata.bind = engine
DBSession = sqlalchemy.orm.sessionmaker(bind=engine)
session = DBSession()

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
    well = i.well
    layout = well.experiment_layout
    mutation_zygosity = 'wt'
    guide_name = 'none'
    
    if i.mutation_summaries:
        mutation_zygosity = i.mutation_summaries[0].zygosity
    if well.well_content.guides: 
        guide_name = well.well_content.guides[0].name
    print(i.sequencing_sample_name, guide_name, mutation_zygosity)
    zygosities.append(mutation_zygosity)
    
print(len(results))

# plots
plotly.offline.plot({
    "data": [go.Histogram(x = zygosities)],
    "layout": go.Layout(title="hello world")
})