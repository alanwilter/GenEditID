# -*- coding: utf-8 -*-
"""
Created on Thu Aug  3 16:13:21 2017

@author: alvare01
"""

import pandas
import sqlalchemy
import logging

import plotly.graph_objs as go
import plotly.offline as py
from plotly import tools

from collections import OrderedDict, defaultdict

from dnascissors.config import cfg
from dnascissors.model import Base
from dnascissors.model import Guide
from dnascissors.model import Target
from dnascissors.model import ExperimentLayout
from dnascissors.model import Project
from dnascissors.model import SequencingLibraryContent
from dnascissors.model import VariantResult
from dnascissors.model import Well
from dnascissors.model import WellContent
from dnascissors.model import MutationSummary

engine = sqlalchemy.create_engine(cfg['DATABASE_URI'])
Base.metadata.bind = engine
DBSession = sqlalchemy.orm.sessionmaker(bind=engine)
dbsession = DBSession()

query = dbsession.query(VariantResult)\
                  .join(VariantResult.sequencing_library_content)\
                  .join(SequencingLibraryContent.well)\
                  .join(SequencingLibraryContent.mutation_summaries)\
                  .join(Well.well_content)\
                  .join(WellContent.guides)\
                  .join(Well.experiment_layout)\
                  .join(ExperimentLayout.project)\
                  .filter(Project.geid == 'GEP00001')\
                  .filter(SequencingLibraryContent.dna_source != 'gDNA')\
                  .filter(VariantResult.allele_fraction > 0.1)\
                  .filter(MutationSummary.has_off_target == False)
results = query.all()
# notes: this query gets only wells with guides, not gDNA, allele fraction > 0.1 and no offtargets

data =[]
guidedict = {} #list of guides with their coordinates
for i in results:
    well = i.sequencing_library_content.well
    variant_caller = i.variant_caller
    plate = well.experiment_layout.geid
    well_position = "{:s}{:02}".format(well.row, well.column)
    guide_name = well.well_content.guides[0].name
    guide_coordinate = well.well_content.guides[0].amplicon_selections[0].guide_location
    mutation_length = 0 if i.indel_length == None else i.indel_length
    if mutation_length ==0:
        mutation = 'SNV'
    elif mutation_length > 0:
        mutation = 'insertion'
    elif mutation_length < 0:
        mutation = 'deletion'
    mutation_coordinate = i.position
    mutation_coordinate_start = mutation_coordinate if mutation_length >=0 else mutation_coordinate - abs(mutation_length)
    mutation_coordinate_end = mutation_coordinate + mutation_length if mutation_length >=0 else mutation_coordinate
    allele_fraction = i.allele_fraction
    data.append([
        variant_caller,
        plate,
        well_position,
        guide_name,
        guide_coordinate % 10000, #to get just the last 4 digits, eg. 12345676 > 5676
        mutation,
        mutation_coordinate_start % 10000,
        mutation_coordinate_end % 10000,
        allele_fraction,
        mutation_length
    ])
    guidedict.update({guide_name:guide_coordinate % 10000})

df = pandas.DataFrame(data)
column_names = ['caller', 'plate', 'well_position', 'guide', 'guide_coordinate',\
 'mutation', 'mutation_coordinate_start', 'mutation_coordinate_end',\
 'allele_fraction', 'indellength']
df = df.rename(columns = dict(enumerate(column_names)))

shapecolor = {'SNV': 'rgba(0, 0, 0, 0.5)', 'insertion': 'rgba(239, 163, 64, 0.5)', 'deletion': 'rgba(112, 161, 239, 0.5)'}
shapetype = {'SNV': 'circle', 'insertion': 'rect', 'deletion': 'rect'}
guidename_log = set()
plotdata_caller = {}
plotdata_xaxisrange = []
plotdata_guides = {}
for a, callerdata in df.groupby(['caller']):
    plotdata_scatter = []
    for i, gDATA in callerdata.groupby(['guide', 'plate', 'well_position']):
        for trow in gDATA.itertuples():
            xaxisrange = [trow.guide_coordinate - 125, trow.guide_coordinate +125]
            xvalues = [trow.mutation_coordinate_start, trow.mutation_coordinate_end]
            yvalues = [trow.allele_fraction] * 2
            hovertext = [', '.join(['well:', '-'.join([trow.plate, trow.well_position]),
                                    'length:', str(trow.indellength),
                                    trow.guide])
                        ] * 2
    #        xvalues = list(range(*[mutation_coordinate_start,trow.mutation_coordinate_end]))
     #       yvalues = trow.allele_fraction * len(xvalues)
            plotdata_scatter.append({
                'x':xvalues,
                'y':yvalues,
                'name':trow.guide,
                'mode':'lines',
                'legendgroup':trow.guide,
                'showlegend':False if trow.guide in guidename_log else True,
                'line':{
                    'color':shapecolor.get(trow.mutation),
                    'width':5
                },
                'text': hovertext
            })
            guidename_log.update([trow.guide])
            plotdata_xaxisrange.append(xaxisrange)
    plotdata_caller.update({trow.caller:plotdata_scatter})


fig = tools.make_subplots(rows = 1, cols = 2, subplot_titles = ['VarDict', 'HaplotypeCaller'])
for i in plotdata_caller['VarDict']:
    fig.append_trace(i, 1, 1)

for j in plotdata_caller['HaplotypeCaller']:
    fig.append_trace(j, 1, 2)

for i in fig.layout:
    if i.find('xaxis') == 0:
        fig.layout[i].update({'range':[min(min(plotdata_xaxisrange)), max(max(plotdata_xaxisrange))], 'showgrid':False})
    elif i.find('yaxis') == 0:
        fig.layout[i].update({'range':[0,1.1]})



#####revise from here, some problem with the shapes (taking only first element of shapes?) and annotations not showing
shapes = []
annots = []
for caller in ['1','2']:
    for guidename,guidecoord in zip(guidedict, guidedict.values()):
        shape={
        'type': 'line',
        'xref':'x'+ caller,
        'yref':'y'+ caller,
        'x0': guidecoord, 'y0': 0, 'x1': guidecoord, 'y1': 1,
        'line': {'color': 'rgb(55, 128, 191)','width': 1, 'dash':'dashdot'},
        }
        annotation={
        'text':guidename,
        'x': guidecoord, 'y':1,
        'xref':'x'+ caller, 'yref':'y'+ caller,
        'textangle':270,
        'font':{'size':9}#,
#        'showarrow':False,
#        'valign':'top',
#        'height':10
        }
        shapes.append(shape)
        annots.append(annotation)
 #   annotations.append(annotation1, annotation2)    

fig.layout.update({
    'shapes': shapes,
    'hovermode': 'closest'})
fig.layout['annotations'].extend(annots)
    
py.plot(fig)


#--------------using shapes instead
#plotdata_scatter = []
#plotdata_shapes = []
#plotdata_xaxisrange = []
#shapecolor = {'SNV': 'black', 'insertion': 'orange', 'deletion': 'steelblue'}
#shapetype = {'SNV': 'circle', 'insertion': 'rect', 'deletion': 'rect'}
#guidename_log = set()
#for i, gDATA in df.groupby(['caller', 'guide', 'plate', 'well_position']):
#    for trow in gDATA.itertuples():
#        xaxisrange = [trow.guide_coordinate - 125, trow.guide_coordinate +125]
#        plotdata_scatter.append({
#            'x':xaxisrange,
#            'y':[0,1],
#            'name':trow.guide,
#            'legendgroup':trow.guide,
#            'showlegend':False if trow.guide in guidename_log else True,
#            'text':'text',
#            'mode':'text'
#        })
#        plotdata_xaxisrange.append(xaxisrange)
#        plotdata_shapes.append({
#            'type': shapetype.get(trow.mutation),
#            'xref': 'x',
#            'yref': 'y',
#            'x0': trow.mutation_coordinate_start,
#            'y0': trow.allele_fraction - 0.2,
#            'x1': trow.mutation_coordinate_end,
#            'y1': trow.allele_fraction,
#            'line': {'width': 0},
#            'fillcolor': shapecolor.get(trow.mutation),
#            'opacity':0.5
#        })
#        guidename_log.update([trow.guide])
#
#plotdata_layout = {
#        'xaxis':{'range':xaxisrange, 'showgrid':False},
#        'yaxis':{'range':[0,1]},
#        'shapes': plotdata_shapes
#        }
#
#
#py.plot(dict(data = plotdata_scatter, layout = plotdata_layout))