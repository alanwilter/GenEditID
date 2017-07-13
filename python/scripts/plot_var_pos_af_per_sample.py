import sqlalchemy
import pandas
import plotly.graph_objs as go
import plotly.offline as pyoff
#import numpy as np

from pandas.core.frame import DataFrame

from dnascissors.config import cfg
from dnascissors.model import Base
from dnascissors.model import SequencingLibraryContent
from dnascissors.model import Well
from dnascissors.model import ExperimentLayout
from dnascissors.model import Project
from dnascissors.model import VariantResult


############################## START OF FUNCTIONS ######################################

def get_ins_shape_coord_string( position, af, length ): # this function is not used
   if length % 2 == 0:
      first_half = length / 2
      first_half = int(first_half)
   else:
      first_half = length / 2
      first_half = int(first_half)
   l1_x = position - first_half
   l1_y = af
   l2_x = position + first_half
   l2_y = af
   ins_coord = 'M ' + str(position) + ' ' + str(0)  +  ' L ' +  ' ' + str(l1_x) + ' ' + ' ' + str(l1_y) + ' L ' + str(l2_x) + ' ' + str(l2_y) + ' Z'
   return ins_coord


def get_guide_shape_coord( guide_pos_list, max_af):
   # vertical line for guide position
   shape_coord_list = []
   for i in guide_pos_list:
      shape_coord = { 'type': 'line', 'x0': i, 'x1': i, 'y0' : 0, 'y1' : max_af,
                      'line': {'color': 'red', 'width': 2, 'dash': 'dashdot' }
                    }
      shape_coord_list.append( shape_coord )
   return shape_coord_list


def get_del_shape_coord( del_start_list, del_len_list, del_af_list ):
   # rectagle box for deletions
   shape_coord_list = []
   for i in range(len(del_start_list)):
      del_start = del_start_list[i]
      del_end = del_start + del_len_list[i]
      del_end -= 1
      y0 = del_af_list[i] - 0.01
      y1 = del_af_list[i]
      shape_coord = {
         'type': 'rect',
         'x0': del_start, 'x1': del_end, 'y0': y0, 'y1': y1,
         'line': { 'width': 0},
         'fillcolor': 'rgba(128, 0, 128, 0.7)',
         'opacity' : 0.5
      }
      shape_coord_list.append( shape_coord )
   return shape_coord_list


def get_ins_shape_coord_triangle(ins_start_list, ins_len_list, ins_af_list): # This function is not used
   # Triange box for deletions
   shape_coord_list = []
   for i in range(len(ins_start_list)):
      shape_string = get_ins_shape_coord_string( position=ins_start_list[i], af=ins_af_list[i], length=ins_len_list[i] )
      shape_coord = {
         'type': 'path',
         'path': shape_string,
         'fillcolor': 'rgba(44, 160, 101, 0.5)',
         'line': {
            'color': 'rgba(44, 160, 101, 0.5)',
            'width': 0
         }
      }
      shape_coord_list.append( shape_coord )
   return shape_coord_list


def get_snv_shape_coord(snv_start_list,  snv_af_list):
   shape_coord_list = []
   for i in range( len(snv_start_list)):
      x0 = snv_start_list[i] - 0.5
      x1 = snv_start_list[i] + 0.5
      y0 = snv_af_list[i] - 0.004
      y1 = snv_af_list[i] + 0.004
      shape_coord = {
         'type': 'circle',
         'x0': x0,  'x1': x1, 'y0': y0, 'y1': y1,
         'fillcolor': 'black',
         'line': { 'color': 'black'}
      }
      shape_coord_list.append(shape_coord)
   return shape_coord_list


def get_ins_shape_coord( ins_start_list, ins_len_list, ins_af_list ):
   # rectagle box for insertions
   shape_coord_list = []
   for i in range(len(ins_start_list)):
      ins_start = ins_start_list[i]
      ins_end = ins_start + ins_len_list[i]
      ins_end -= 1
      y0 = ins_af_list[i] - 0.01
      y1 = ins_af_list[i]
      shape_coord = {
         'type': 'rect',
         'x0': ins_start, 'x1': ins_end, 'y0': y0, 'y1': y1,
         'line': { 'width': 0},
         'fillcolor': 'blue',
         'opacity' : 0.5
      }
      shape_coord_list.append( shape_coord )
   return shape_coord_list


def get_last_3_digits( pos_list):
   pos_list_trim = []
   for i in pos_list:
      pos_list_trim.append( i % 1000)
   return pos_list_trim

############################## END OF FUNCTIONS #########################################

# sqlalchemy query

engine = sqlalchemy.create_engine(cfg['DATABASE_URI'])
Base.metadata.bind = engine
DBSession = sqlalchemy.orm.sessionmaker(bind=engine)
session = DBSession()

results = session.query(SequencingLibraryContent)\
                  .join(SequencingLibraryContent.variant_results)\
                  .join(SequencingLibraryContent.well)\
                  .join(Well.well_content)\
                  .join(Well.experiment_layout)\
                  .join(ExperimentLayout.project)\
                  .filter(Project.geid == 'GEP00001')\
                  .all()

mutations = []
for sequencing_library_content in results:
   well = sequencing_library_content.well
   layout = well.experiment_layout
   variant_caller = 'VarDict'
   allele_fraction_threshold = 0

   for variant in sequencing_library_content.variant_results:
      # select INDEL variants only from VarDict caller with allele_fraction above 0.1
      if variant.variant_caller == variant_caller and variant.variant_type == 'INDEL' and variant.allele_fraction > allele_fraction_threshold:
         if len(well.well_content.guides) == 1:
            # find which cut sites match the variant
            guide = well.well_content.guides[0]
            matched_amplicon_selection = None
            for amplicon_selection in guide.amplicon_selections:
               if variant.chromosome.endswith(amplicon_selection.amplicon.chromosome) and \
                               variant.position >= amplicon_selection.amplicon.start and \
                               variant.position <= amplicon_selection.amplicon.end and amplicon_selection.is_on_target:
                  matched_amplicon_selection = amplicon_selection
            if matched_amplicon_selection:
               # select mutation only if within the amplicon range
               row = [variant.variant_caller, \
                      variant.variant_type, \
                      variant.position, \
                      variant.allele_fraction, \
                      variant.indel_length, \
                      sequencing_library_content.sequencing_sample_name, \
                      well.well_content.clone.name, \
                      guide.name, \
                      matched_amplicon_selection.guide_location, \
                      matched_amplicon_selection.is_on_target, \
                      '{}_chr{}_{}'.format(matched_amplicon_selection.amplicon.genome.assembly, matched_amplicon_selection.amplicon.chromosome, matched_amplicon_selection.amplicon.start), \
                      matched_amplicon_selection.amplicon.start, \
                      matched_amplicon_selection.amplicon.end]
               mutations.append(row)

#for m in mutations:
#   print(m)

df = DataFrame(mutations, columns=['variant_caller', 'variant_type', 'position', \
                                   'allele_fraction', 'indel_length', 'sequencing_sample_name',\
                                   'clone.name', 'name', 'guide_location', 'is_on_target',\
                                   '', 'start', 'end'])

#print( df)
col_names = df.columns.values
col_names[10] = 'amplicon'

#dfgroup = df.groupby(['amplicon'])
dfgroup = df.groupby(['sequencing_sample_name'])

#print(dfgroup)
#print(col_names )

for sample_name, grouped_data in dfgroup:
   #print( sample_name )
   #print( grouped_data )

   # amplicon start and end positions
   # positions are too bing to show on the plot
   # therefore take last three digits
   # Threfore take 100 bp up and 100 bp down stream of amplicon for plotting
   #
   amp_start = get_last_3_digits(grouped_data['start'])[0] -100
   amp_end = get_last_3_digits(grouped_data['end'])[0] + 100

   # Guide positions
   max_af = 1
   guide_pos = list(set(grouped_data['guide_location']) )
   clone_name = list(set(grouped_data['clone.name']))
   #print( guide_pos )

   # get varinat start, length and af
   var_type_data = list( grouped_data['variant_type'] )
   var_position = list( grouped_data['position'] )

   #print( var_position)
   var_af = list(grouped_data['allele_fraction'])
   # convert alle fraction to percentage??
   #var_af = list(np.array( var_af) * 100 )
   #print( type(var_af[0]) )
   indel_length = list(grouped_data['indel_length'])
   ins_start = []
   ins_length = []
   ins_af = []
   del_start = []
   del_length = []
   del_af = []
   snv_pos =[]
   snv_af = []

   for j in range( len(var_type_data)):
      if var_type_data[j] == 'INDEL':
         # sort into insertions and deletions
         if indel_length[j] < 0:
            del_length.append(abs(indel_length[j]))
            del_start.append( var_position[j] )
            del_af.append( var_af[j] )

         if indel_length[j] > 0:
            ins_length.append(indel_length[j])
            ins_start.append( var_position[j] )
            ins_af.append( var_af[j] )
      if var_type_data[j] == 'SNV':
         snv_pos.append( var_position[j])
         snv_af.append( var_af[j] )

   trace0 = go.Scatter(
      x=[amp_start, amp_end],
      y=[0, 1],
      mode='text'
   )
   data = [trace0]

   guide_pos = get_last_3_digits(guide_pos)
   ins_start = get_last_3_digits(ins_start)
   del_start = get_last_3_digits(del_start)
   snv_pos = get_last_3_digits(snv_pos)

   guide_lo_list = get_guide_shape_coord(guide_pos_list=guide_pos, max_af=max_af)
   ins_lo_list = get_ins_shape_coord(ins_start_list=ins_start, ins_len_list=ins_length, ins_af_list=ins_af)
   del_lo_list = get_del_shape_coord(del_start_list=del_start, del_len_list=del_length, del_af_list=del_af)
   snv_lo_list = get_snv_shape_coord( snv_start_list=snv_pos, snv_af_list=snv_af)

   shapes_layout_list = guide_lo_list + ins_lo_list + del_lo_list + snv_lo_list

   layout = {
      'xaxis': {'range': [amp_start, amp_end], 'showgrid': False, 'title' : 'Position', 'visible' : True },
      'yaxis': {'range': [0, max_af], 'showgrid': False, 'title' : 'Allele Fraction', 'visible' : True},
      'shapes': shapes_layout_list,
      'title' : sample_name
   }

   fig_file_name = sample_name + '.html'
   fig = {
      'data': data,
      'layout': layout
   }
   pyoff.plot(fig, filename=fig_file_name)

   #break