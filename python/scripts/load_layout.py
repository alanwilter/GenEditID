import sqlalchemy

from geneditid.config import cfg
from geneditid.model import Base
from geneditid.loader import ProjectConfigLoader

import log as logger


def main():
    import argparse
    import os
    parser = argparse.ArgumentParser()
    parser.add_argument("--layout", dest="file_layout", action="store", help="The file layout e.g. '20170118_GEP00001.xlsx'", required=True)
    parser.add_argument("--geid", dest="geid", action="store", help="Project GEP ID e.g. 'GEP00001'", required=True)
    parser.add_argument("--clean", dest="clean_db", action="store_true", default=False, help="Clean database before loading?")
    options = parser.parse_args()

    log = logger.get_custom_logger(os.path.join(cfg['PROJECTS_FOLDER'], 'load_layout.log'))

    engine = sqlalchemy.create_engine(cfg['DATABASE_URI'])
    Base.metadata.bind = engine
    DBSession = sqlalchemy.orm.sessionmaker(bind=engine)
    session = DBSession()

    loader = ProjectConfigLoader(session, options.geid, options.file_layout)

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
