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
from dnascissors.model import ExperimentLayout
from dnascissors.model import Well
from dnascissors.model import Project
from dnascissors.model import SequencingLibraryContent
from dnascissors.model import VariantResult

#-------- sqlalchemy query
engine = sqlalchemy.create_engine(cfg['DATABASE_URI'])
Base.metadata.bind = engine
DBSession = sqlalchemy.orm.sessionmaker(bind=engine)
session = DBSession()

results = session.query(VariantResult)\
                 .join(VariantResult.sequencing_library_content)\
                 .join(SequencingLibraryContent.well)\
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
marker_symbol = ['circle', 'triangle-up', 'cross', 'hash']
data = []
nloop = -1
for caller, gdata in df.groupby(['caller']):
    nloop += 1
    for guidename, gdata2 in gdata.groupby(['guides']):
        gdata2_byvar = gdata2.groupby(['indellengths']).size()
        gdata2_byvar_percent = gdata2_byvar*100 / gdata2_byvar.sum()
        htext = []
        print(gdata2_byvar_percent)
        print(len(gdata2))
        print(gdata2_byvar_percent.tolist())
        for length_gdbp in range(len(gdata2_byvar_percent)):
            hovertext = [caller] #, guidename, '%=' + str(round(gdata2_byvar_percent.tolist()[length_gdbp]))] In case we want to add other things to the hover data
            hovertext = ' '.join(hovertext)
            htext.append(hovertext)
        #for index, rows in gdata2_byvar_percent():
          #  print(index,rows)
            #htext = [caller, guidename, '%=' + str(round(gdata2_byvar_percent, 2))]
           # htext = ', '.join(htext) # the hover text needs to be a single string, it can't be a list
           # hovertext.append(htext)
           # print(hovertext)
        data.append(
            go.Scatter(
                 x = gdata2_byvar_percent.index.tolist(),
                 y = gdata2_byvar_percent.tolist(),
                 name = guidename,
                 mode = 'markers',
                 marker = {'symbol': marker_symbol[nloop]},
                 text = htext
            )
        )


layout = go.Layout(
    title = 'Indel lengths',
    yaxis = {'title': '% of submitted samples per guide'}    
)

# plot
pyoff.plot(dict(data = data, layout = layout))

