# -*- coding: utf-8 -*-
"""
Created on Wed Jun 28 14:27:18 2017

@author: alvare01

plot description: 
Query database to get guides and zygosities ('wt', 'homo', 'dmut'...).
Calculation of percentages of zygosities ('wt', 'homo', 'dmut'...) on a per-guide basis
Bar plot

"""

import sqlalchemy
import plotly.graph_objs as go
import plotly.offline as pyoff

from dnascissors.model import MutationSummary
from dnascissors.model import Base
from dnascissors.config import cfg
from dnascissors.model import SequencingLibraryContent
from dnascissors.model import Well
from dnascissors.model import ExperimentLayout
from dnascissors.model import Project
from pandas import DataFrame
from pandas import groupby

#-------- sqlalchemy query
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

guide = []
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
    guide.append(guide_name)
    
print(len(results))

#-------- calculate percentages per guide in a pandas dataframe
# convert 'results' to pandas dataframe and group by 'guides'
df = DataFrame({'guides': guide, 'zygosities': zygosities})
dfgroup = df.groupby(['guides'])

#for guidename, grouped_data in dfgroup:
#    print(guidename)
#    grouped_data_byzygosity = grouped_data.groupby(['zygosities']).size()
#    grouped_data_byzygosity_percent = grouped_data_byzygosity*100 / grouped_data_byzygosity.sum()
#    print(grouped_data_byzygosity_percent)
#    
#for guidename in dfgroup:
#    grouped_data_byzygosity = grouped_data.groupby(['zygosities']).size()
#    grouped_data_byzygosity_percent = grouped_data_byzygosity*100 / grouped_data_byzygosity.sum()

# calculate percentages of zygosities and create bar plot 'data' dictionary
data = []
for guidename, grouped_data in dfgroup:
    grouped_data_byzygosity = grouped_data.groupby(['zygosities']).size()
    grouped_data_byzygosity_percent = grouped_data_byzygosity*100 / grouped_data_byzygosity.sum()
    data.append(
        go.Bar(
             x = grouped_data_byzygosity_percent.index.tolist(),
             y = grouped_data_byzygosity_percent.tolist(),
             name = guidename
             )
    )
    

# order x axis values
vals = ['wt', 'homo', 'smut', 'dmut', 'iffy']
layout = go.Layout(
    title = 'Zygosities',
    xaxis = {'categoryorder': 'array', 'categoryarray': vals},
    yaxis = {'title': '% of submitted samples per guide'}    
)

# plot
pyoff.plot(dict(data = data, layout = layout))
