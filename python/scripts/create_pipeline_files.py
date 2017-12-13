import sqlalchemy

from dnascissors.config import cfg
from dnascissors.model import *
from dnascissors.loader import CellGrowthLoader
from shutil import copyfile

import log as logger

def convert_strand(s):
    if s == 'forward': return '+'
    if s == 'reverse': return '-'
    return '*'

def create_amplicon_file(session, project, seq_dict, driver_file):
    copyfile(seq_dict, driver_file)
    out = open(driver_file, "a")
    try:
        selectionQuery = session.query(AmpliconSelection)\
                            .join(AmpliconSelection.amplicon)\
                            .join(AmpliconSelection.guide)\
                            .join(Guide.target)\
                            .join(Target.project)\
                            .filter(Project.geid == project)\
                            .order_by(Amplicon.chromosome.asc(), Amplicon.start.asc())
        
        for s in selectionQuery.all():
            a = s.amplicon
            out.write("{}\t{}\t{}\t{}\t{}\n".format(a.chromosome, a.start, a.end, convert_strand(s.guide_strand),\
                                                    "{}_{}_{}".format(a.genome.assembly, a.chromosome, a.start)))
    finally:
        out.close()

def create_target_file(session, project, seq_dict, driver_file):
    copyfile(seq_dict, driver_file)
    out = open(driver_file, "a")
    try:
        targetQuery = session.query(Target)\
                             .join(Target.project)\
                             .filter(Project.geid == project)\
                             .order_by(Target.chromosome.asc(), Target.start.asc())
        
        for t in targetQuery.all():
            out.write("{}\t{}\t{}\t{}\t{}\n".format(t.chromosome, t.start, t.end, convert_strand(t.strand),\
                                                    "{}_{}_{}".format(t.genome.assembly, t.chromosome, t.start)))
    finally:
        out.close()

def main():
    import argparse
    import os

    parser = argparse.ArgumentParser()
    parser.add_argument("--project", dest="project", action="store", help="The project id.", required=True)
    parser.add_argument("--seq-dict", dest="dict", action="store", help="The reference sequence dictionary.", required=True)
    options = parser.parse_args()

    log = logger.get_custom_logger(os.path.join(os.path.dirname(__file__), 'create_pipeline_files.log'))

    engine = sqlalchemy.create_engine(cfg['DATABASE_URI'])
    Base.metadata.bind = engine
    DBSession = sqlalchemy.orm.sessionmaker(bind=engine)
    session = DBSession()

    try:
        create_amplicon_file(session, options.project, options.dict, "amplicons.txt")
        create_target_file(session, options.project, options.dict, "targets.txt")
        
        session.commit()
    except Exception as e:
        log.exception(e)
        session.rollback()
    finally:
        session.close()


if __name__ == '__main__':
    main()
