import sqlalchemy

from dnascissors.config import cfg
from dnascissors.model import Base
from dnascissors.loader import VariantLoader

import log as logger


def main():
    import argparse
    import os

    parser = argparse.ArgumentParser()
    parser.add_argument("--file", dest="file", action="store", help="The variant raw results excel file e.g. 'data/20170127_GEP00001/SLX-13775.vardict.variants.xlsx'.", required=True)
    parser.add_argument("--project_geid", dest="project_geid", action="store", help="The project geid e.g. 'GEP00001'", required=True)
    parser.add_argument("--caller", dest="variant_caller", action="store", help="The type of variant caller used.", choices=['VarDict', 'HaplotypeCaller', 'FreeBayes'], required=True)
    parser.add_argument("--clean", dest="clean_db", action="store_true", default=False, help="Clean database before loading?")
    options = parser.parse_args()

    log = logger.get_custom_logger(os.path.join(os.path.dirname(__file__), 'load_variant_results.log'))

    engine = sqlalchemy.create_engine(cfg['DATABASE_URI'])
    Base.metadata.bind = engine
    DBSession = sqlalchemy.orm.sessionmaker(bind=engine)
    session = DBSession()
    loader = VariantLoader(session, options.project_geid, options.file, options.variant_caller)

    try:
        if options.clean_db:
            loader.clean()
        loader.load()
        session.commit()
    except Exception as e:
        log.exception(e)
        session.rollback()
    finally:
        session.close()


if __name__ == '__main__':
    main()
