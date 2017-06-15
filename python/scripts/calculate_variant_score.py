import os
import sqlalchemy
from dnascissors.config import cfg
from dnascissors.model import Base
from dnascissors.model import VariantResult
from dnascissors.model import SequencingLibraryContent
from dnascissors.model import Well
from dnascissors.model import ExperimentLayout
from dnascissors.model import Project
import log as logger


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
    Base.metadata.bind = engine
    DBSession = sqlalchemy.orm.sessionmaker(bind=engine)
    dbsession = DBSession()
    results = dbsession.query(VariantResult)\
                       .join(VariantResult.sequencing_library_content)\
                       .join(SequencingLibraryContent.well)\
                       .join(Well.well_content)\
                       .join(Well.experiment_layout)\
                       .join(ExperimentLayout.project)\
                       .filter(Project.geid == 'GEP00001')\
                       .filter(VariantResult.allele_fraction > 0.15)\
                       .filter(VariantResult.allele_fraction < 0.90)\
                       .filter(sqlalchemy.or_(VariantResult.sift.like('tolerated%'), VariantResult.sift == None))\
                       .all()

    for result in results:
        log.info('----------')
        log.info('variant: {} {}'.format(result.chromosome, result.position))
        sequencing_library_content = result.sequencing_library_content
        well = sequencing_library_content.well
        # check amplicon and allele are on the same chromosome
        if not result.amplicon.split('_')[1] == result.chromosome:
            raise ValueError('Amplicon on {:s} and variant on {:s}'.format(result.amplicon.split('_')[1], result.chromosome))
        if len(well.well_content.guides) == 1:
            # find which cut sites match the variant
            guide = well.well_content.guides[0]
            matched_amplicon_selection = None
            for amplicon_selection in guide.amplicon_selections:
                if result.chromosome.endswith(amplicon_selection.amplicon.chromosome) and result.position >= amplicon_selection.amplicon.start and result.position <= amplicon_selection.amplicon.end:
                    matched_amplicon_selection = amplicon_selection
                    # calculate the distance between cut site and variant position
                    distance = abs(int(result.position) - int(amplicon_selection.guide_location))
                    # calculate score based on distance
                    if distance == 0:
                        score = 1
                    elif distance <= 5 and distance > 0:
                        score = 0.75
                    elif distance < 10 and distance > 5:
                        score = 0.50
                    else:
                        score = 0.25
            if matched_amplicon_selection:
                # store result in database
                log.info('guide: {}, onTarget: {}, cutSite: {}, distance: {}, score: {}'.format(guide.name, matched_amplicon_selection.is_on_target, matched_amplicon_selection.guide_location, distance, score))
            else:
                log.info('guide: {}, not match found for this variant (position outside start/end of all amplicons), no score calculated'.format(guide.name))
        elif len(well.well_content.guides) > 1:
            raise Exception('More than one associated guide, {} found. Cannot calculate the score.'.format(len(well.well_content.guides)))
        else:
            log.info('No guide found, no score calculated')


if __name__ == '__main__':
    main()
