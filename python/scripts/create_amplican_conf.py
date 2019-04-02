import sqlalchemy
from dnascissors.config import cfg
from dnascissors.model import Base
import log as logger
import amplifind


def create_files(session, refgenome, project):
    amplican_file = "amplican_config.csv"

    with open(amplican_file, "w") as output:
        output.write("ID,Barcode,guideRNA,Forward_Primer,Reverse_Primer,Direction,Amplicon\n")
        for sample in amplifind.get_samples(session, refgenome, project):
            output.write('{},{},{},{},{},{},{}\n'.format(sample['name'],
                                                         sample['barcode'],
                                                         sample['guide'],
                                                         sample['fprimer'],
                                                         sample['rprimer'],
                                                         sample['direction'],
                                                         sample['amplicon']))


def main():
    import argparse
    import os

    parser = argparse.ArgumentParser()
    parser.add_argument("--project", dest="project", action="store", help="The project id.", required=True)
    parser.add_argument("--genome", dest="refgenome", action="store", help="The reference genome fasta file e.g. 'hsa.GRCh38_hs38d1.fa'", required=True)
    options = parser.parse_args()

    log = logger.get_custom_logger(os.path.join(os.path.dirname(__file__), 'create_amplican_conf.log'))

    engine = sqlalchemy.create_engine(cfg['DATABASE_URI'])
    Base.metadata.bind = engine
    DBSession = sqlalchemy.orm.sessionmaker(bind=engine)
    session = DBSession()

    try:
        create_files(session, options.refgenome, options.project)
        session.commit()
    except Exception as e:
        log.exception(e)
        session.rollback()
    finally:
        session.close()


if __name__ == '__main__':
    main()
