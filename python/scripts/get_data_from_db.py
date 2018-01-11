import sqlalchemy
import csv

from dnascissors.config import cfg
from dnascissors.model import Base
from dnascissors.model import SequencingLibraryContent
from dnascissors.model import Well
from dnascissors.model import ExperimentLayout
from dnascissors.model import Project

# sqlalchemy database connection
engine = sqlalchemy.create_engine(cfg['DATABASE_URI'])
Base.metadata.bind = engine
DBSession = sqlalchemy.orm.sessionmaker(bind=engine)
session = DBSession()

# sqlalchemy query
results = session.query(SequencingLibraryContent)\
                 .outerjoin(SequencingLibraryContent.mutation_summaries)\
                 .outerjoin(SequencingLibraryContent.variant_results)\
                 .join(SequencingLibraryContent.well)\
                 .join(Well.well_content)\
                 .join(Well.experiment_layout)\
                 .join(ExperimentLayout.project)\
                 .filter(Project.geid == 'GEP00001')\
                 .all()

# Growth data: Unique ID | hours | confluence_percentage
# Protein data: Unique ID | intensity_channel_700 | intensity_channel_800
# Variant results: Unique ID | variant_caller | variant_type | chromosome | position | ref | alt | allele_fraction | depth | amplicon | indel_length
# Mutation summaries: Unique ID | zygosity | consequence

data = []
for sequencing_library_content in results:
    well = sequencing_library_content.well
    row = {'well_id': well.id,
           'well_location': '{}_{}{}'.format(well.experiment_layout.geid, well.row, well.column),
           'sample_name': sequencing_library_content.sequencing_sample_name,
           'abundance_700': well.abundances[0].intensity_channel_700,
           'abundance_800': well.abundances[0].intensity_channel_800,
           'growths': [(g.hours, g.confluence_percentage) for g in well.growths],
           'variant_results': [(v.variant_caller, v.variant_type, v.chromosome, v.position, v.ref, v.alt, v.allele_fraction, v.depth, v.amplicon, v.indel_length) for v in sequencing_library_content.variant_results],
           'mutation_summaries': [(m.zygosity, m.consequence, m.variant_caller_presence, m.score) for m in sequencing_library_content.mutation_summaries]
          }
    data.append(row)

with open('GEP00001_abundances.csv', 'w') as csv_abundances,\
     open('GEP00001_growths.csv', 'w') as csv_growths,\
     open('GEP00001_mutation_summaries.csv', 'w') as csv_mutation_summaries,\
     open('GEP00001_variant_results.csv', 'w') as csv_variant_results:
    writer_csv_abundances = csv.DictWriter(csv_abundances, fieldnames=['well_id', 'well_location', 'sample_name', 'abundance_700', 'abundance_800'], extrasaction='ignore', delimiter='\t')
    writer_csv_abundances.writeheader()
    writer_csv_growths = csv.writer(csv_growths, delimiter='\t')
    writer_csv_growths.writerow(['well_id', 'well_location', 'sample_name'] + ['growth_{}'.format(g[0]) for g in data[0]['growths']])
    writer_csv_mutation_summaries = csv.writer(csv_mutation_summaries, delimiter='\t')
    writer_csv_mutation_summaries.writerow(['well_id', 'well_location', 'sample_name', 'zygosity', 'merged_consequence', 'score'])
    writer_variant_results = csv.writer(csv_variant_results, delimiter='\t')
    writer_variant_results.writerow(['well_id', 'well_location', 'sample_name', 'variant_caller', 'variant_type', 'chromosome', 'position', 'ref', 'alt', 'allele_fraction', 'depth', 'amplicon', 'indel_length'])
    for d in data:
        writer_csv_abundances.writerow(d)
        writer_csv_growths.writerow([d['well_id'], d['well_location'], d['sample_name']] + [g[1] for g in d['growths']])
        if len(d['mutation_summaries']) > 0:
            writer_csv_mutation_summaries.writerow([d['well_id'], d['well_location'], d['sample_name'], d['mutation_summaries'][0][0], d['mutation_summaries'][0][1], d['mutation_summaries'][0][3]])
        for variant_result in d['variant_results']:
            writer_variant_results.writerow([d['well_id'], d['well_location'], d['sample_name']] + [v for v in variant_result])
