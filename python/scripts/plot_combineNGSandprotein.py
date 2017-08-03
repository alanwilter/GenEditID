# -*- coding: utf-8 -*-
"""
Created on Thu Aug  3 11:05:52 2017

@author: alvare01
"""
import pandas
import sqlalchemy
import logging

import plotly.graph_objs as go
import plotly.offline as py
from plotly import tools

from dnascissors.config import cfg
from dnascissors.model import Base
from dnascissors.model import ExperimentLayout
from dnascissors.model import Project
from dnascissors.model import Well
from dnascissors.model import WellContent

color_scale = [[0.0, 'rgb(185, 190, 193)'],  # e.g. this color is applied to the first 10% of values
               [0.3, 'rgb(200, 223, 247)'],
               [0.5, 'rgb(237, 119, 90)'],
               [0.6, 'rgb(237, 195, 90)'],
               [0.7, 'rgb((237, 229, 90)'],
               [0.8, 'rgb((181, 221, 24)'],
               [1, 'rgb(39, 132, 1)']]
               
               
engine = sqlalchemy.create_engine(cfg['DATABASE_URI'])
Base.metadata.bind = engine
DBSession = sqlalchemy.orm.sessionmaker(bind=engine)
dbsession = DBSession()

query = dbsession.query(Well)\
              .join(Well.experiment_layout)\
              .join(ExperimentLayout.project)\
              .join(Well.well_content)\
              .join(WellContent.guides)\
              .filter(Project.geid == 'GEP00001')
#note I'm getting only the results with guides, ignoring the controls, etc.            
results = query.all()
if not results:
    return None

# get all guides in the project
guidenames = []
for well in results:
    guidenames.append(well.well_content.guides[0].name)

guidenames = set(guidenames)

# build a list of plot traces
data = []
for well in results:
    zygosity = 'notsequenced'
    consequence = 'notsequenced'
    score = 0
    well_abundance_ratio = 0
    
    if well.abundances:
        well_abundance_ratio = well.abundances[0].ratio_800_700
    if well.well_content:
        well_content_type = 'sample-well'
        zygosity = 'wt'
        guide_name = well.well_content.guides[0].name
    if well.sequencing_library_contents:
        # TODO need a loop here to make sure dna_source is fixed cells
        # cannot trust order of items in list returned by sqlalchemy
        # sequencing_library_contents[0].dna_source corresponds to fixed cells,
        # and sequencing_library_contents[1].dna_source to gDNA
        consequence = 'wt'
        if well.sequencing_library_contents[0].mutation_summaries:
            score = well.sequencing_library_contents[0].mutation_summaries[0].score
            zygosity = well.sequencing_library_contents[0].mutation_summaries[0].zygosity
            consequence = well.sequencing_library_contents[0].mutation_summaries[0].consequence
    data_loop = ([
        "{:s}{:02}".format(well.row, well.column),
        well.experiment_layout.geid,
        consequence,
        str(round(well_abundance_ratio, 3)),
        guide_name,
        score,
        zygosity
    ])
    htext = [', '.join([str(i) for i in data_loop])]
    data.append(data_loop + htext)


df = pandas.DataFrame(data)
column_names = ['wellposition', 'plate', 'consequence', 'protein', 'guide', 'score', 'zygosity', 'hovertext']
df = df.rename(columns=dict(enumerate(column_names)))
# group by plate layout
dfgroup = df.groupby(['guide'])

dataplot = []
for i, grouped_data in dfgroup:
    scatter = {
        'mode': 'markers',
        'x': grouped_data.consequence.tolist(),
        'y': grouped_data.protein.tolist(),
        'name': i,
        'legendgroup': i,
        'marker': {
            'color':grouped_data.score.tolist(),  # assign a color based on score
            'colorscale':color_scale,
            'cmin':7000,  # min value of the colorscale
            'cmax':9000,  # max value of the colorscale
            'colorbar':{  # side color bar custom size, fraction of plot size (otherwise it takes all plot height)
                'lenmode': "fraction",
                'len': 0.4
            }
         },
        'text': grouped_data.hovertext.tolist()
    }
    dataplot.append(scatter)

categories = sorted(set(df.consequence))
print(categories)
layout = go.Layout(
    title='Combined NGS and protein',
    xaxis={'categoryorder': 'array', 'categoryarray': categories},
    yaxis={'title': 'Relative protein abundance'},
    hovermode = 'closest'
)

py.plot(dict(data = dataplot, layout = layout))