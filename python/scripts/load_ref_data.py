import os
import sqlalchemy

from dnascissors.config import cfg
from dnascissors.model import Base
from dnascissors.loader import RefLoader

import log as logger


def main():

    log = logger.get_custom_logger(os.path.join(cfg['PROJECTS_FOLDER'], 'load_ref_data.log'))

    engine = sqlalchemy.create_engine(cfg['DATABASE_URI'])
    Base.metadata.bind = engine
    DBSession = sqlalchemy.orm.sessionmaker(bind=engine)
    session = DBSession()

    loader = RefLoader(session)

    try:
        loader.load_genomes()
        loader.load_celllines()
        session.commit()
    except Exception as e:
        log.exception(e)
        session.rollback()
    finally:
        session.close()


if __name__ == '__main__':
    main()
