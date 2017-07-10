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
import pandas
import plotly.graph_objs as go
import plotly.offline as pyoff

from dnascissors.config import cfg
from dnascissors.model import Base
from dnascissors.model import SequencingLibraryContent
from dnascissors.model import Well
from dnascissors.model import ExperimentLayout
from dnascissors.model import Project
from dnascissors.model import VariantResult

#-------- sqlalchemy query
engine = sqlalchemy.create_engine(cfg['DATABASE_URI'])
Base.metadata.bind = engine
DBSession = sqlalchemy.orm.sessionmaker(bind=engine)
session = DBSession()

results = session.query(VariantResult)\
                 .join(VariantResult.sequencing_library_content)\
                 .join(SequencingLibraryContent.well)\
                 .join(Well.well_content)\
                 .join(Well.experiment_layout)\
                 .join(ExperimentLayout.project)\
                 .filter(Project.geid == 'GEP00001')\
                 .filter(SequencingLibraryContent.dna_source != 'gDNA')\
                 .filter(VariantResult.allele_fraction > 0.1)\
                 .all()
                 
guide = []
allindellengths = []
typecaller = []

for i in results:
    well = i.sequencing_library_content.well
    layout = well.experiment_layout
    indellength = 0
    guide_name = 'none'
    
    if well.well_content.guides: 
        guide_name = well.well_content.guides[0].name
    print(i.sequencing_library_content.sequencing_sample_name, guide_name, i.indel_length, i.variant_caller)
    guide.append(guide_name)    
    allindellengths.append(i.indel_length)
    typecaller.append(i.variant_caller)

print(len(results))


#-------- calculate percentages per guide in a pandas dataframe
# convert 'results' to pandas dataframe and group by 'guides'
df = pandas.DataFrame({'caller': typecaller, 'guides': guide, 'indellengths': allindellengths})

# calculate percentages of indel lengths and create bar plot 'data' dictionary
data = []
for caller, gdata in df.groupby(['caller']):
    for guidename, gdata2 in gdata.groupby(['guides']):
        gdata2_byvar = gdata.groupby(['indellengths']).size()
        gdata2_byvar_percent = gdata2_byvar*100 / gdata2_byvar.sum()
        data.append(
            go.Bar(
            x = gdata2_byvar_percent.index.tolist(),
            y = gdata2_byvar_percent.tolist(),
            legendgroup = caller,
            name = guidename
            )
        )


layout = go.Layout(
    title = 'Indel lengths',
    yaxis = {'title': '% of submitted samples per guide'}    
)

# plot
pyoff.plot(dict(data = data, layout = layout))

