import sqlalchemy

from dnascissors.config import cfg
from dnascissors.model import *
from shutil import copyfile

import log as logger


def convert_strand(s):
    if s == 'forward':
        return '+'
    if s == 'reverse':
        return '-'
    return '*'


def create_files(session, project, seq_dict):
    amp_file = "amplicons.txt"
    target_file = "targets.txt"
    
    copyfile(seq_dict, amp_file)
    copyfile(seq_dict, target_file)
    
    with open(amp_file, "a") as amp_out, open(target_file, "a") as tar_out:
        ampliconQuery = session.query(Amplicon)\
                            .join(AmpliconSelection.amplicon)\
                            .join(AmpliconSelection.guide)\
                            .join(Guide.target)\
                            .join(Target.project)\
                            .filter(Project.geid == project)\
                            .order_by(Amplicon.chromosome.asc(), Amplicon.start.asc())
        
        # Should do this as another query with IN [amplicon.in_(Primer.amplicons], but it doesn't seem to work.
        primers = session.query(Primer).all()

        for a in ampliconQuery.all():
            
            #print("Amplicon {} {} {}-{}".format(a.id, a.chromosome, a.start, a.end))
                  
            amp_out.write("{}\t{}\t{}\t+\t{}\n".format(a.chromosome, a.start, a.end,
                                                        "{}_{}_{}".format(a.genome.assembly, a.chromosome, a.start)))
            
            aprimers = [ p for p in primers if a in p.amplicons ]
            
            tstart = a.start + 1
            tend = a.end - 1
                
            for p in aprimers:
                #print("Has primer {} {}-{} {} {}".format(p.id, p.start, p.end, p.strand, p.sequence))
                
                if p.strand == 'forward':
                    tstart = p.start + 1
                if p.strand == 'reverse':
                    tend = p.end - 1
            
            tar_out.write("{}\t{}\t{}\t+\t{}\n".format(a.chromosome, tstart, tend,
                                                        "{}_{}_{}".format(a.genome.assembly, a.chromosome, a.start)))


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
        create_files(session, options.project, options.dict)
        session.commit()
    except Exception as e:
        log.exception(e)
        session.rollback()
    finally:
        session.close()


if __name__ == '__main__':
    main()
