# -*- coding: utf-8 -*-
"""
Created on Thu Aug  3 16:13:21 2017

@author: alvare01
"""

import sqlalchemy

import plotly.offline as py
from plotly import tools

from dnascissors.config import cfg
from dnascissors.model import Base
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

# get the caller and sequencing_sample_names
results_callers = list(set(i.variant_caller for i in query.all()))
results_samples = list(set(i.sequencing_library_content.sequencing_sample_name for i in query.all()))
# initialise the plot dictionaries
plotdict = {}
for caller in results_callers:
    plotcaller = {}
    for sam in results_samples:
        plotcaller.update({sam:[]})
    plotdict.update({caller:plotcaller})

# populate the plot dictionaries on the fly from the database
data =[]
guidedict = {} #list of guides with their coordinates
guidename_log = set()
plotdata_xaxisrange = []
shapecolor = {'SNV': 'rgba(0, 0, 0, 0.5)', 'insertion': 'rgba(239, 163, 64, 0.5)', 'deletion': 'rgba(112, 161, 239, 0.5)'}
shapetype = {'SNV': 'circle', 'insertion': 'rect', 'deletion': 'rect'}
for i in results:
    well = i.sequencing_library_content.well
    variant_caller = i.variant_caller
    seqsample = i.sequencing_library_content.sequencing_sample_name
    plate = well.experiment_layout.geid
    well_position = "{:s}{:02}".format(well.row, well.column)
    guide_name = well.well_content.guides[0].name
    guide_coordinate = well.well_content.guides[0].amplicon_selections[0].guide_location % 10000 #to get just the last 4 digits, eg. 12345676 > 5676
    xaxisrange = [guide_coordinate - 125, guide_coordinate +125]
    mutation_length = 0 if i.indel_length == None else i.indel_length
    if mutation_length ==0:
        mutation = 'SNV'
    elif mutation_length > 0:
        mutation = 'insertion'
    elif mutation_length < 0:
        mutation = 'deletion'
    mutation_coordinate = i.position
    mutation_coordinate_start = mutation_coordinate if mutation_length >=0 else mutation_coordinate - abs(mutation_length)
    mutation_coordinate_start =mutation_coordinate_start % 10000
    mutation_coordinate_end = mutation_coordinate + mutation_length if mutation_length >=0 else mutation_coordinate
    mutation_coordinate_end = mutation_coordinate_end % 10000
    allele_fraction = i.allele_fraction
    hovertext = [', '.join([variant_caller,'well: ' + '-'.join([plate, well_position]),
                            'length: '+ str(mutation_length), guide_name])] * 2
    guidedict.update({guide_name:guide_coordinate})
    plotdict[variant_caller][seqsample].append({
        'x':[mutation_coordinate_start, mutation_coordinate_end],
        'y':[allele_fraction]*2,
        'name':guide_name,
        'mode':'lines',
        'legendgroup':guide_name,
        'showlegend':False if guide_name in guidename_log else True,
        'line':{
            'color':shapecolor.get(mutation),
            'width':5
        },
        'text': hovertext  
    }) 
    guidename_log.update([guide_name])
    plotdata_xaxisrange.append(xaxisrange)
    

fig = tools.make_subplots(rows = 1, cols = 2, subplot_titles = results_callers)
col_index = 1
for d_caller in plotdict:
    for d_sample in plotdict[d_caller]:
        for mut in plotdict[d_caller][d_sample]:
            fig.append_trace(mut, 1, col_index)
    col_index += 1

for i in fig.layout:
    if i.find('xaxis') == 0:
        fig.layout[i].update({'range':[min(min(plotdata_xaxisrange)), max(max(plotdata_xaxisrange))], 'showgrid':False})
    elif i.find('yaxis') == 0:
        fig.layout[i].update({'range':[0,1.1]})

shapes,annots = [],[]
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