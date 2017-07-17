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

from plotly import tools
from math import ceil, floor
from itertools import cycle

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
results = session.query(Well)\
                 .join(Well.well_content)\
                 .join(Well.experiment_layout)\
                 .join(ExperimentLayout.project)\
                 .outerjoin(SequencingLibraryContent.variant_results)\
                 .outerjoin(SequencingLibraryContent.sequencing_sample_name)\
                 .filter(Project.geid == 'GEP00001')\
                 .filter(SequencingLibraryContent.dna_source != 'gDNA')\
                 .all()
                
#                query(SequencingLibraryContent)\
#                 .outerjoin(SequencingLibraryContent.mutation_summaries)\
#                 .join(SequencingLibraryContent.well)\
#                 .join(Well.well_content)\
#                 .join(Well.experiment_layout)\
#                 .join(ExperimentLayout.project)\
#                 .filter(Project.geid == 'GEP00001')\
#                 .filter(SequencingLibraryContent.dna_source != 'gDNA')\
#                 .all()
    
    
    
df = []
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
    df.append(rowzygo)

# convert to pandas dataframe
df = pandas.DataFrame(df)
# column_names = ['sample', 'plate', 'row', 'column', 'guide', 'protein_abundance', 'zygosity']
dfgroup = df.groupby([1])



















# plot
data = []
for i, grouped_data in dfgroup:
    xnull = None # apparently I have to have values in all x and y coordinates so the plot axes show all coordinates (1:12 and A:H)
    ynull = None
    scatter =   go.Scatter(
                    x=grouped_data[2] if grouped_data[2] else x = xnyll,
                    y=grouped_data[3] if grouped_data[3],
                    name=i,
                    mode = 'markers',
                    marker = dict(
                        size = '20',
                        color = df[df[4] == 'STAT3.2'][5],
                        colorscale='Viridis',
                        showscale=True
                        )
                    )
    data.append(scatter)


#trace0 = go.Scatter(
#    x = df[df[4] == 'STAT3.2'][2],
#    y = df[df[4] == 'STAT3.2'][3],
#    mode = 'markers',
#    marker = dict(
#        size = '20',
#        color = df[df[4] == 'STAT3.2'][5],
#        colorscale='Viridis',
#        showscale=True
#        )
#    )

#To remove a trace's colorbar use the showscale41 attribute.
#
#in your case, you could remove all but the first trace's colorbar and make that colorbar span [0, 50] using cmin12 and cmax3.
    
#layout = go.Layout(
#    autosize=False,
#    width=400,
#    height=230,
#    margin=go.Margin(
#        l=15,
#        r=15,
#        b=15,
#        t=15
#    )
#)


# create plot layout in a grid of two columns and n rows
 # calculate number of subplots (number of plates)
numberofplates = len(set(df[1]))
numberofplotrows = ceil(numberofplates / 2)

# create figure
fig = tools.make_subplots(rows=numberofplotrows, cols=2)
row_counter = 1
for i,j in zip(range(0, numberofplates), data):
    subplot_row    = floor(row_counter)  # the roughest way to change row every two loops? Oh dear.
    row_counter    = row_counter + 0.5
    subplot_column = 1 if 1&i == 0 else 2
    print(row_counter, subplot_row, subplot_column)
    fig.append_trace(j, subplot_row, subplot_column)

# layout construction. In the scatter subplots, each plot has an independent axis 
  # (xaxis, xaxis2... xaxisN, and the same for yaxis). Each layout for each subplot
  # must be built independently.

for i in fig.layout:
    fig.layout[i].update({'type': 'category', 'autoscale':'False'})
    if i.find('xaxis') == 0:
        fig.layout[i].update({'categoryorder': 'array', 'categoryarray': list(range(1,13))})
    elif i.find('yaxis') == 0:
        fig.layout[i].update({'categoryorder': 'array', 'categoryarray': ['A', 'B', 'C', 'D', 'E', 'F', 'G']})

pyoff.plot(fig)




layout = go.Layout(
    xaxis = dict(range = list(range(1,13))),
    yaxis = dict(range = [['A', 'B', 'C', 'D', 'E', 'F', 'G']])
)

#the idea is to loop through the keys in the layout, and if 'xaxis' is found add {'range': list(range(1,13))},
# and if 'yaxis' is found add {'range': ['A', 'B', 'C', 'D', 'E', 'F', 'G']} using the .update method.
# example a.update({'range':['A', 'B', 'C']})
for i in fig.layout:
    for n in i.keys():
        n.find('xaxis')
 

   
fig['layout'].update(
    xaxis = dict(range = list(range(1,13))),
    yaxis = dict(range = [['A', 'B', 'C', 'D', 'E', 'F', 'G']])
    )

#fig.append_trace(data[0], 1, 1)
#fig.append_trace(data[1], 1, 2)
#fig.append_trace(data[2], 2, 1)
#fig.append_trace(data[3], 2, 2)
pyoff.plot(fig)







plotlayout = []
for i in range(numberofplates):
    layout = go.Layout(        
        xaxis = {'domain': [0, subplot_xposition]},
        yaxis = {'domain': [0, subplot_xposition]}
        )
    plotlayout.append(layout)

fig = tools.make_subplots(rows=1, cols=2)

fig.append_trace(trace1, 1, 1)
fig.append_trace(trace2, 1, 2)
    
    
layout = []



layout = go.Layout(
    xaxis=dict(
        domain=[0, 0.45 * 1&i]
    ),
    yaxis=dict(
        domain=[0, 0.45]

