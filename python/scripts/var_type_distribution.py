import sqlalchemy
import pandas
import plotly.graph_objs as go
import plotly.offline as pyoff
# import numpy as np

from pandas.core.frame import DataFrame

from dnascissors.config import cfg
from dnascissors.model import Base
from dnascissors.model import SequencingLibraryContent
from dnascissors.model import Well
from dnascissors.model import ExperimentLayout
from dnascissors.model import Project
from dnascissors.model import VariantResult

# sqlalchemy query

engine = sqlalchemy.create_engine(cfg['DATABASE_URI'])
Base.metadata.bind = engine
DBSession = sqlalchemy.orm.sessionmaker(bind=engine)
session = DBSession()

results = session.query(SequencingLibraryContent) \
   .join(SequencingLibraryContent.variant_results) \
   .join(SequencingLibraryContent.well) \
   .join(Well.well_content) \
   .join(Well.experiment_layout) \
   .join(ExperimentLayout.project) \
   .filter(Project.geid == 'GEP00001') \
   .all()

mutations = []
for sequencing_library_content in results:
   well = sequencing_library_content.well
   layout = well.experiment_layout
   #variant_caller = 'VarDict'
   allele_fraction_threshold = 0

   for variant in sequencing_library_content.variant_results:
      # select INDEL variants only from VarDict caller with allele_fraction above 0.1
      if  variant.variant_type == 'INDEL' and variant.allele_fraction > allele_fraction_threshold:
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
                      variant.consequence, \
                      variant.position, \
                      variant.allele_fraction, \
                      variant.indel_length, \
                      sequencing_library_content.sequencing_sample_name, \
                      well.well_content.clone.name, \
                      guide.name, \
                      matched_amplicon_selection.guide_location, \
                      matched_amplicon_selection.is_on_target, \
                      '{}_chr{}_{}'.format(matched_amplicon_selection.amplicon.genome.assembly,
                                           matched_amplicon_selection.amplicon.chromosome,
                                           matched_amplicon_selection.amplicon.start), \
                      matched_amplicon_selection.amplicon.start, \
                      matched_amplicon_selection.amplicon.end]
               mutations.append(row)


df = DataFrame(mutations, columns=['variant_caller', 'variant_type', 'consequence', 'position', \
                                   'allele_fraction', 'indel_length', 'sequencing_sample_name', \
                                   'clone.name', 'name', 'guide_location', 'is_on_target', \
                                   '', 'start', 'end'])

col_names = df.columns.values
col_names[10] = 'amplicon'

# dfgroup = df.groupby(['amplicon'])
dfgroup = df.groupby(['variant_caller'])

bar_data = []
for variant_caller, grouped_data in dfgroup:
   var_type_data = list(grouped_data['variant_type'])
   indel_length = list(grouped_data['indel_length'])
   snv_counts, ins_counts, del_counts = 0,0,0


   for j in range(len(var_type_data)):
      if var_type_data[j] == 'SNV':
         snv_counts += 1
      if var_type_data[j] == 'INDEL':
         if indel_length[j] > 0:
            ins_counts += 1
         if indel_length[j] < 0:
            del_counts += 1

   total_var_counts = snv_counts + ins_counts + del_counts

   snv_per = ( snv_counts / total_var_counts) * 100
   ind_per = ( ins_counts / total_var_counts) * 100
   del_per = ( del_counts / total_var_counts) * 100

   bar_data.append(
      go.Bar(
         x=['SNVs', 'Insertions', 'Deletions'],
         y=[snv_per, ind_per, del_per],
         name=variant_caller
      )
   )

bar_layout = {
   'title' : 'Variant Type Distribution',
   'barmode' : 'group',
   'yaxis': { 'title': '% Variants','range': [0,100] }
}

fig = {
   'data' : bar_data,
   'layout' : bar_layout
}

pyoff.plot(fig)
