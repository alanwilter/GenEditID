# -*- coding: utf-8 -*-
"""
Created on Mon Jul 10 18:32:07 2017

@author: alvare01

Plot description:
Query database to get variant callers, guides and alleles (e.g. C/CATG).
Calculation of percentages of alleles on a per-guide and per-allele basis
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
from collections import OrderedDict

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
alleles = []
typecaller = []
allele_index5 = [] #this will give the position of '/' in the allele, to then sort the dataframe
allele_index3 = []

for i in results:
    well = i.sequencing_library_content.well
    layout = well.experiment_layout
    indellength = 0
    guide_name = 'none'
    
    if well.well_content.guides: 
        guide_name = well.well_content.guides[0].name
    print(i.sequencing_library_content.sequencing_sample_name, guide_name, i.alleles, i.variant_caller)
    guide.append(guide_name)    
    alleles.append(i.alleles)
    typecaller.append(i.variant_caller)
    allele_split = i.alleles.split('/')
    allele_index5.append(len(allele_split[0]))
    allele_index3.append(len(allele_split[1]))

print(len(results))


#-------- calculate percentages per guide in a pandas dataframe
# convert 'results' to pandas dataframe and group by 'guides'
df = pandas.DataFrame({'caller': typecaller, 'guides': guide, 'alleles': alleles, 'index5': allele_index5, 'index3': allele_index3})
df = df.sort_values(by = ['index5', 'index3']) # calculate percentages of indel lengths and create bar plot 'data' dictionary

data = []
for caller, gdata in df.groupby(['caller']):
    for guidename, gdata2 in gdata.groupby(['guides']):
        gdata2_byvar = gdata2.groupby(['alleles']).size()
        gdata2_byvar_percent = gdata2_byvar*100 / gdata2_byvar.sum()
        data.append(
            go.Bar(
            x = gdata2_byvar_percent.tolist(),
            y = gdata2_byvar_percent.index.tolist(),
            legendgroup = caller,
            name = guidename,
            orientation = 'h'
            )
        )

# order the y-axis
 # get the ordered alleles in the y-axis
allele_list_sorted = list(OrderedDict.fromkeys(df.alleles))

layout = go.Layout(
    title = 'Alleles',
    yaxis = {'title': '% of alleles in submitted samples per guide',
             'categoryorder': 'array',
             'categoryarray': allele_list_sorted} ,
    margin = {'l': 800}
)

# plot
pyoff.plot(dict(data = data, layout = layout))

