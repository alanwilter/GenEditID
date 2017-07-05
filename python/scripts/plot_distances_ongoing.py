# -*- coding: utf-8 -*-
"""
Created on Wed Jul  5 14:54:56 2017

@author: alvare01

plot description: 
Query database to get guides and indel lengths.
Calculation of percentages of indel lengths ('wt', 'homo', 'dmut'...) on a per-guide and per-allele basis
Bar plot

"""

import sqlalchemy
import plotly.graph_objs as go
import plotly.offline as pyoff

# from dnascissors.model import MutationSummary
from dnascissors.model import Base
from dnascissors.model import SequencingLibraryContent
from dnascissors.model import Well
from dnascissors.model import ExperimentLayout
from dnascissors.model import Project
from dnascissors.model import VariantResult

from dnascissors.config import cfg

from pandas import DataFrame
from pandas import groupby

#-------- sqlalchemy query
engine = sqlalchemy.create_engine(cfg['DATABASE_URI'])
Base.metadata.bind = engine
DBSession = sqlalchemy.orm.sessionmaker(bind=engine)
session = DBSession()

results = session.query(SequencingLibraryContent)\
                 .outerjoin(SequencingLibraryContent.variant_results)\
                 .join(SequencingLibraryContent.well)\
                 .join(Well.well_content)\
                 .join(Well.experiment_layout)\
                 .join(ExperimentLayout.project)\
                 .filter(Project.geid == 'GEP00001')\
                 .filter(SequencingLibraryContent.dna_source != 'gDNA')\
                 .all()
                 
guide = []
allindellengths = []
for i in results:
    well = i.well
    layout = well.experiment_layout
    indellength = 0
    guide_name = 'none'
    
    if i.variant_results:
        indellength = i.variant_results[0].indel_length
    if well.well_content.guides: 
        guide_name = well.well_content.guides[0].name
    print(i.sequencing_sample_name, guide_name, indellength)
    allindellengths.append(indellength)
    guide.append(guide_name)

print(len(results))

#-------- calculate percentages per guide in a pandas dataframe
# convert 'results' to pandas dataframe and group by 'guides'
df = DataFrame({'guides': guide, 'indellengths': allindellengths}) ##PROBLEM HERE! I expect a higher number of lenghts (this is per allele), and don't know why I'm getting NaN's
dfgroup = df.groupby(['guides'])

# calculate percentages of indel lengths and create bar plot 'data' dictionary
data = []
for guidename, grouped_data in dfgroup:
    grouped_data_byvar = grouped_data.groupby(['indellengths']).size()
    grouped_data_byvar_percent = grouped_data_byvar*100 / grouped_data_byvar.sum()
    data.append(
        go.Bar(
             x = grouped_data_byvar_percent.index.tolist(),
             y = grouped_data_byvar_percent.tolist(),
             name = guidename
             )
    )
    

# order x axis values
layout = go.Layout(
    title = 'Indel lengths',
    yaxis = {'title': '% of submitted samples per guide'}    
)

# plot
pyoff.plot(dict(data = data, layout = layout))
