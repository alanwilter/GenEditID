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
import pandas
import plotly.graph_objs as go
import plotly.offline as pyoff

from dnascissors.config import cfg
from dnascissors.model import Base
from dnascissors.model import ExperimentLayout
from dnascissors.model import Project
from dnascissors.model import SequencingLibraryContent
from dnascissors.model import VariantResult
from dnascissors.model import Well

# from plotly import tools. Possible use to make subplots

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
typemutation = []
typevariant = []
typecaller = []

for i in results:
    well = i.sequencing_library_content.well
    layout = well.experiment_layout
    typemut = 'wt'
    guide_name = 'none'
    
    if well.well_content.guides:
        guide_name = well.well_content.guides[0].name
    print(guide_name, i.consequence, i.variant_type, i.variant_caller)
    
    guide.append(guide_name)
    typemutation.append(i.consequence)
    typevariant.append(i.variant_type)
    typecaller.append(i.variant_caller)

print(len(results))


#-------- calculate percentages per guide in a pandas dataframe
# convert 'results' to pandas dataframe and group by 'variants'
df = pandas.DataFrame({'variants': typevariant, 'caller': typecaller, 'guides': guide, 'mutation': typemutation})
dfgroup_typevariant = df.groupby(['variants'])

# create independent datasets (one for indels, one for snvs) for plotting
dfgroup_typevariant_INDEL = dfgroup_typevariant.get_group('INDEL')  #[i for i in dfgroup_typevariant][0][1]
dfgroup_typevariant_SNV = dfgroup_typevariant.get_group('SNV')      #[i for i in dfgroup_typevariant][1][1]

# calculate percentages of mutation types and create bar plot 'data' dictionary
# This will produce a plot with grouped legends. Annoyingly there no feature at date 20170710 to add 
#  titles to the legend groups.
#  It's been suggested to use layout annotations as a workaround: https://github.com/plotly/plotly.js/issues/689
#  I assume that the top legend group is the first variant caller.
marker_symbol = ['circle', 'triangle-up', 'cross', 'hash']

data_INDELs = []
nloop = -1
for caller, gdata in dfgroup_typevariant_INDEL.groupby(['caller']):
    nloop +=1
    for guidename, gdata2 in gdata.groupby(['guides']):
#        print(caller)
#        print(guidename)
#        print(gdata2)
        gdata2_byvar = gdata2.groupby(['mutation']).size()
        gdata2_byvar_percent = gdata2_byvar*100 / gdata2_byvar.sum()
        data_INDELs.append(
            go.Scatter(
                 x = gdata2_byvar_percent.index.tolist(),
                 y = gdata2_byvar_percent.tolist(),
                 name = guidename,
                 mode = 'markers',
                 marker = {'symbol':  marker_symbol[nloop]},
                 text = caller
                 )
        ) 

data_SNVs = []
nloop = -1
#marker_color = ['rgb(0,0,0)', 'rgb(230,159.0)', 'rgb(86,180,233)',
#                'rgb(0.158,115)', 'rgb(240,228,66)', 'rgb(0,114,178)',
#                'rgb(213,94,0)', 'rgb(204, 121,167)']
for caller, gdata in dfgroup_typevariant_SNV.groupby(['caller']):
    nloop +=1
    for guidename, gdata2 in gdata.groupby(['guides']):
#        print(caller)
#        print(guidename)
#        print(gdata2)
        gdata2_byvar = gdata2.groupby(['mutation']).size()
        gdata2_byvar_percent = gdata2_byvar*100 / gdata2_byvar.sum()
        data_SNVs.append(
            go.Scatter(
                 x = gdata2_byvar_percent.index.tolist(),
                 y = gdata2_byvar_percent.tolist(),
                 name = guidename,
                 mode = 'markers',
                 marker = {'symbol':  marker_symbol[nloop]},
                 text = caller
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