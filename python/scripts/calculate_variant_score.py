import os
import sqlalchemy
from dnascissors.config import cfg
import log as logger
import pandas

"""
filtering step:
---------------
  1. filter on allele frequency
    - keep alleles with frequency > 0.15 and < 0.90 ??? or above 90 ???
    - most interesting alleles have a frequency between 0.4 and 0.6

  2. filter out on SIFT with value starting with 'tolerated'

scoring step:
-------------
  1. score based on distance from cut site
  distance calculated using Position and guide_location from amplicon_selection_table
  transform into percentage?

  2. add on/off target information
  filtering out off-target?

  x. score based on indel length, score decreases when length increases — talk to Ruben
  not yet in the scoring calculation

  x. score down: inframe_deletion
  not yet in the scoring calculation

"""


def main():
    log = logger.get_custom_logger(os.path.join(os.path.dirname(__file__), 'calculate_variant_score.log'))
    engine = sqlalchemy.create_engine(cfg['DATABASE_URI'])
    query = """
   with guide_names as (
      select
         guide_well_content_association.well_content_id,
         string_agg(guide.name, ', ') as guide_names
      from
         guide_well_content_association
         inner join guide on guide_well_content_association.guide_id = guide.id
      group by guide_well_content_association.well_content_id
   ),

   cut_sites as (
      select
         guide_well_content_association.well_content_id,
         string_agg(concat('chr', amplicon.chromosome, '_', amplicon_selection.guide_location, '_', amplicon_selection.is_on_target), ', ') as cut_sites
      from
         guide_well_content_association
         inner join guide on guide_well_content_association.guide_id = guide.id
         inner join amplicon_selection on amplicon_selection.guide_id = guide.id
         inner join amplicon on amplicon.id = amplicon_selection.amplicon_id
      group by guide_well_content_association.well_content_id
   )

  select
      experiment_layout.geid as experiment_layout_geid,
      concat(well.row, well.column) as well,
      cell_line.name as cell_line_name,
      clone.name as clone_name,
      guide_names.guide_names,
      cut_sites.cut_sites,
      well_content.is_control,
      well_content.content_type,
      sequencing_library_content.dna_source,
      sequencing_library_content.sequencing_sample_name,
      sequencing_library_content.sequencing_barcode,
      sequencing_library.slxid,
      variant_raw_result.*
   from
      well
      inner join experiment_layout on experiment_layout.id = well.experiment_layout_id
      inner join project on project.id = experiment_layout.project_id
      left join well_content on well_content.id = well.well_content_id
      left join guide_names on guide_names.well_content_id = well_content.id
      left join cut_sites on cut_sites.well_content_id = well_content.id
      left join clone on clone.id = well_content.clone_id
      left join cell_line on cell_line.id = clone.cell_line_id
      left join sequencing_library_content on sequencing_library_content.well_id = well.id
      inner join sequencing_library on sequencing_library.id = sequencing_library_content.sequencing_library_id
      inner join variant_raw_result on variant_raw_result.sequencing_library_content_id = sequencing_library_content.id
   where
      project.geid = 'GEP00001'
      and variant_raw_result.allele_fraction > 0.15
      and variant_raw_result.allele_fraction < 0.90
      and not variant_raw_result.sift = 'tolerated%%'

    """
    with engine.connect() as con:
        df = pandas.read_sql(query, con)
        for row in df.itertuples():
            # check amplicon and variant are on the same chromosome
            variant_chr = row.chromosome
            amplicon_genome, amplicon_chr, amplicon_start = row.amplicon.split('_')
            if not amplicon_chr == variant_chr:
                raise ValueError('Amplicon on {:s} and variant on {:s}'.format(amplicon_chr, variant_chr))
            if row.cut_sites:
                # create list of cut sites
                cut_sites = row.cut_sites.split(', ')
                # find which cut site match variant_chr
                the_cut_site_chr, the_cut_site_position, the_cut_site_on_target = None, None, None
                for cut_site in cut_sites:
                    cut_site_chr, cut_site_position, cut_site_on_target = cut_site.split('_')
                    if cut_site_chr == variant_chr:
                        the_cut_site_chr, the_cut_site_position, the_cut_site_on_target = cut_site_chr, cut_site_position, cut_site_on_target
                # report on/off target information
                # report cut site position
                # calculate the distance between cut site and variant positions
                distance = abs(int(row.position) - int(the_cut_site_position))
                # calculate score based on distance
                if distance == 0:
                    score = 1
                elif distance <= 5 and distance > 0:
                    score = 0.75
                elif distance < 10 and distance > 5:
                    score = 0.50
                else:
                    score = 0.25
                log.debug('onTarget: {:s}, cutSite: {:s}, distance: {:d}, score: {:f}'.format(the_cut_site_on_target, the_cut_site_position, distance, score))


if __name__ == '__main__':
    main()
