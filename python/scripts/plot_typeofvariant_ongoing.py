# -*- coding: utf-8 -*-
"""
Created on Wed Jul  5 15:34:57 2017

@author: alvare01

plot description: 
Query database to get guides and type of mutation ('wt', 'frameshift', 'inframe deletion'...).
Calculation of percentages of types of mutation on a per-guide and per-allele basis
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
# from plotly import tools. Possible use to make subplots

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
                 #.filter(SequencingLibraryContent.variant_results[0].allele_fraction > 0.1)\  this filter needs to be added
guide = []
typemutation = []
typevariant = []
for i in results:
    well = i.well
    layout = well.experiment_layout
    typemut = 'wt'
    guide_name = 'none'
    
    if i.variant_results:
        typemut = i.variant_results[0].consequence
        typevar = i.variant_results[0].variant_type
    if well.well_content.guides: 
        guide_name = well.well_content.guides[0].name
    print(i.sequencing_sample_name, guide_name, typemut, typevar)
    typemutation.append(typemut)
    typevariant.append(typevar)
    guide.append(guide_name)

print(len(results))

#-------- calculate percentages per guide in a pandas dataframe
# convert 'results' to pandas dataframe and group by 'variants'
df = DataFrame({'variants': typevariant, 'guides': guide, 'mutation': typemutation}) ##PROBLEM HERE! I expect a higher number of lenghts (this is per allele), and don't know why I'm getting NaN's
dfgroup_typevariant = df.groupby(['variants'])

# create independent datasets (one for indels, one for snvs) for plotting
dfgroup_typevariant_INDEL = [i for i in dfgroup_typevariant][0][1]
dfgroup_typevariant_SNV = [i for i in dfgroup_typevariant][1][1]

# calculate percentages of mutation types and create bar plot 'data' dictionary
data_INDELs = []
for guidename, grouped_data in dfgroup_typevariant_INDEL.groupby(['guides']):
    grouped_data_byvar = grouped_data.groupby(['mutation']).size()
    grouped_data_byvar_percent = grouped_data_byvar*100 / grouped_data_byvar.sum()
    data_INDELs.append(
        go.Bar(
             x = grouped_data_byvar_percent.index.tolist(),
             y = grouped_data_byvar_percent.tolist(),
             name = guidename
             )
    )

data_SNVs = []
for guidename, grouped_data in dfgroup_typevariant_SNV.groupby(['guides']):
    grouped_data_byvar = grouped_data.groupby(['mutation']).size()
    grouped_data_byvar_percent = grouped_data_byvar*100 / grouped_data_byvar.sum()
    data_SNVs.append(
        go.Bar(
             x = grouped_data_byvar_percent.index.tolist(),
             y = grouped_data_byvar_percent.tolist(),
             name = guidename
             )
    )

 
# layouts
layout_INDELs = go.Layout(
    title = 'Type of mutation (INDELS)',
    yaxis = {'title': '% of submitted samples per guide'}    
)

layout_SNVs = go.Layout(
    title = 'Type of mutation (SNVs)',
    yaxis = {'title': '% of submitted samples per guide'}    
)

# plot
pyoff.plot(dict(data = data_INDELs, layout = layout_INDELs))
pyoff.plot(dict(data = data_SNVs, layout = layout_SNVs))
# not sure why, but when the script is sourced, both plots are the SNV plot


# the ideal situation would be to plot both datasets side by side
 # pyoff.plot(dict(data = data_INDELs + data_SNVs, layout = layout_INDELs)) #say we use a different layout for each dataset
 # however the guides appear duplicated. As far as I have seen, you can't have 
 # the same legent for two different plots in plotly (ironically you can have it in if you use ggplotly(R:ggplot2)... )
