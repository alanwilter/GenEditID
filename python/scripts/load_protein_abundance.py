import sqlalchemy

from dnascissors.config import cfg
from dnascissors.model import Base
from dnascissors.loader import ProteinAbundanceLoader

import log as logger


def main():
    import argparse
    import os

    parser = argparse.ArgumentParser()
    parser.add_argument("--plateid", dest="plateid", action="store", help="The plate ID e.g. 'GEP00001_01'", required=True)
    parser.add_argument("--file", dest="file", action="store", help="The InCell Western file e.g. 'data/20170127_GEP00001/GEP00001_01_ICW.csv'.", required=True)
    options = parser.parse_args()

    log = logger.get_custom_logger(os.path.join(cfg['PROJECTS_FOLDER'], 'load_protein_abundance.log'))

    engine = sqlalchemy.create_engine(cfg['DATABASE_URI'])
    Base.metadata.bind = engine
    DBSession = sqlalchemy.orm.sessionmaker(bind=engine)
    session = DBSession()
    loader = ProteinAbundanceLoader(session, options.file, options.plateid)

    try:
        loader.load()
        session.commit()
    except Exception as e:
        log.exception(e)
        session.rollback()
    finally:
        session.close()


if __name__ == '__main__':
    main()
