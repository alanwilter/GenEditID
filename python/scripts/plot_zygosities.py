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
import pandas
import plotly.graph_objs as go
import plotly.offline as py

from dnascissors.config import cfg
from dnascissors.model import Base
from dnascissors.model import ExperimentLayout
from dnascissors.model import Project
from dnascissors.model import SequencingLibraryContent
from dnascissors.model import Well


# sqlalchemy database connection
engine = sqlalchemy.create_engine(cfg['DATABASE_URI'])
Base.metadata.bind = engine
DBSession = sqlalchemy.orm.sessionmaker(bind=engine)
session = DBSession()

# sqlalchemy query
#we need to filter by gDNA, because there is a sample in project1 (GE-P6B4-G) 
 #that is gDNA only (it was sent for sequencing only as gDNA, without a 'fixed cells' counterpart)
results = session.query(Well)\
        .join(Well.sequencing_library_contents)\
        .join(Well.experiment_layout)\
        .join(ExperimentLayout.project)\
        .filter(Project.geid == 'GEP00001')\
        .filter(SequencingLibraryContent.dna_source != 'gDNA')\
        .all()


guides = []
zygosities = []
for well in results:
    #well = i.well
    layout = well.experiment_layout
    mutation_zygosity = 'wt'
    guide_name = 'none'
    if well.sequencing_library_contents[0].mutation_summaries:
        mutation_zygosity = well.sequencing_library_contents[0].mutation_summaries[0].zygosity
    if well.well_content.guides:
        guide_name = well.well_content.guides[0].name
    #print(well.sequencing_library_contents[0].sequencing_sample_name, guide_name, mutation_zygosity)
    zygosities.append(mutation_zygosity)
    guides.append(guide_name)

# convert 'results' to pandas dataframe and group by 'guides'
df = pandas.DataFrame({'guides': guides, 'zygosities': zygosities})
dfgroup = df.groupby(['guides'])

# calculate percentages of zygosities and create bar plot 'data' dictionary
plots = []
for guide_name, grouped_data in dfgroup:
    grouped_data_byzygosity = grouped_data.groupby(['zygosities']).size()
    grouped_data_byzygosity_percent = grouped_data_byzygosity*100 / grouped_data_byzygosity.sum()
    plots.append(
        go.Bar(
             x=grouped_data_byzygosity_percent.index.tolist(),
             y=grouped_data_byzygosity_percent.tolist(),
             name=guide_name
             )
    )

# order x axis values
categories = ['wt', 'homo', 'smut', 'dmut', 'iffy']


print(categories)
layout = go.Layout(
    title='Zygosities',
    xaxis={'categoryorder': 'array', 'categoryarray': categories},
    yaxis={'title': '% of submitted samples per guide'}
)

# plot
py.plot(dict(data=plots, layout=layout))
