# -*- coding: utf-8 -*-
"""
Created on Tue Jul 11 12:32:07 2017

@author: alvare01

plot description: 
Query database to get sample names, plate, well locations, guides, scores, protein value, slope value and zygosity per clone
Create a scatter plot with the shape of a 96-well plate, mapping the scores to 
    a colour gradient and showing hover data: guide, zygosity, protein ratio and score

"""

import sqlalchemy
import pandas
import plotly.graph_objs as go
import plotly.offline as pyoff

from plotly import tools
from math import ceil, floor

from dnascissors.config import cfg
from dnascissors.model import Base
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
                 .join(Well.experiment_layout)\
                 .join(ExperimentLayout.project)\
                 .filter(Project.geid == 'GEP00001')\
                 .all()

# data collection
data = []
for i, well in enumerate(results, start=1):
    zygosity = 'empty-well'
    guide_name = 'no-guide'
    sequencing_sample_name = 'not-sequenced'
    dna_source = 'not-sequenced'
    well_content_type = 'empty-well'
    score = 0
    row = [well.experiment_layout.project.geid,
           well.row,
           well.column,
           well.abundances[0].ratio_800_700
           ]
    if well.well_content:
        well_content_type = 'sample-well'
        zygosity = 'wt'
        if well.well_content.guides:
            guide_name = well.well_content.guides[0].name
    if well.sequencing_library_contents:
        sequencing_sample_name = well.sequencing_library_contents[0].sequencing_sample_name
        if well.sequencing_library_contents[0].mutation_summaries:
            dna_source = well.sequencing_library_contents[0].dna_source
            score = well.sequencing_library_contents[0].mutation_summaries[0].score
            # sequencing_library_contents[0].dna_source corresponds to fixed cells,
            # and sequencing_library_contents[1].dna_source to gDNA
            zygosity = well.sequencing_library_contents[0].mutation_summaries[0].zygosity
    data.append([
                well.experiment_layout.project.geid,
                well.experiment_layout.geid,
                well.row,
                well.column,
                guide_name,
                well.abundances[0].ratio_800_700,
                zygosity,
                score,
                well_content_type
                ])

# convert to pandas dataframe and add column names
df = pandas.DataFrame(data)
column_names = ['project', 'plate', 'row', 'column', 'guide', 'protein_abundance', 'zygosity', 'score', 'well_content_type']
df = df.rename(columns = dict(enumerate(column_names)))

# group by plate layout
dfgroup = df.groupby(['plate'])

## plot

 # custom colorscale for the score color bar. The scale goes from 0 to 1
 # (0% to 100% of values in the plot)
color_scale=[[0.0, 'rgb(239, 239, 239)'], # e.g. this color is applied to the first 10% of values
            [0.3, 'rgb(200, 223, 247)'],
            [0.5, 'rgb(237, 119, 90)'],
            [0.6, 'rgb(237, 195, 90)'],
            [0.7, 'rgb((237, 229, 90)'],
            [0.8, 'rgb((181, 221, 24)'],
            [1, 'rgb(39, 132, 1)']]
            
 # construct the data dictionary            
dataplot = []
nloop = 0
for i, grouped_data in dfgroup:
    nloop += 1
    # create hover text to add to the scatter data structure
    hovertext = []
    for index, rows in grouped_data.iterrows():
        htext = [rows.guide, rows.zygosity, 'protein=' + str(round(rows.protein_abundance, 2)), 'score=' + str(rows.score)]
        htext = ', '.join(htext) # the hover text needs to be a single string, it can't be a list
        hovertext.append(htext)
    # create scatter plot data structure
    scatter = go.Scatter(
        x = grouped_data['column'],
        y = grouped_data['row'],
        name = i,        # i is the plate
        mode = 'markers',
        marker = dict(
            size = '20',                   # dot size
            color = grouped_data['score'], # assign a color based on score
            showscale = True if nloop == 1 else False,
            colorscale = color_scale,
            cmin = 7000, # min value of the colorscale
            cmax = 9000, # max value of the colorscale
            colorbar = { # side color bar custom size, fraction of plot size (otherwise it takes all plot height)
                'lenmode':"fraction",
                'len':0.4                
            }
        ),
        # hover text
        text = hovertext
    )
    dataplot.append(scatter)

# create plot layout in a grid of two columns and n rows
 # calculate number of subplots (number of plates)
numberofplates = len(set(df['plate']))
numberofplotrows = ceil(numberofplates / 2)
plotheight = numberofplotrows*330 # calculation of total plot height, see comment further down

# create figure
subplottitles = [i for i,j in dfgroup]
fig = tools.make_subplots(rows = numberofplotrows, cols = 2, subplot_titles = subplottitles)
row_counter = 1
for i,j in zip(range(0, numberofplates), dataplot):
    subplot_row    = floor(row_counter)  # the roughest way to change row every two loops? Oh dear.
    row_counter    = row_counter + 0.5
    subplot_column = 1 if 1&i == 0 else 2
    print(row_counter, subplot_row, subplot_column)
    fig.append_trace(j, subplot_row, subplot_column)

# update layout construction. In the scatter subplots, each plot has an independent axis 
  # (xaxis, xaxis2... xaxisN, and the same for yaxis). Each layout for each subplot
  # must be built independently.

for i in fig.layout:
    if i.find('xaxis') == 0:
        fig.layout[i].update({'categoryorder': 'array', 'categoryarray': list(range(1,13))}) # this is to show all xaxis values in the plot, in this order
        fig.layout[i].update({'type': 'category'})
    elif i.find('yaxis') == 0:
        fig.layout[i].update({'categoryorder': 'array', 'categoryarray': ['H', 'G', 'F', 'E', 'D', 'C', 'B', 'A']}) # this is to show all yaxis values in the plot, in this order
        fig.layout[i].update({'type': 'category'})
        fig.layout[i].update({'type': 'category'})

# with this layout.update below, width and height are set in the plot.Perhaps you can set them directly on the plotting area on the web page
    # hovermode = closest shows the values for the hover point, otherwise by default ('compare' mode) you only see one coordinate
fig.layout.update(dict(autosize = False, width = 850, height = plotheight, hovermode = 'closest'))

# plot
pyoff.plot(fig, filename = 'plot96.html')