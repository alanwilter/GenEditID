import os

import geneditid.log as logger
from geneditid.config import cfg
from geneditid.connect import dbsession
from geneditid.loader import RefLoader


def main():

    log = logger.get_custom_logger(os.path.join(cfg['PROJECTS_FOLDER'], 'load_ref_data.log'))

    try:
        loader = RefLoader(dbsession)
        loader.load_genomes()
        loader.load_celllines()
        dbsession.commit()
    except Exception as e:
        log.exception(e)
        dbsession.rollback()
    finally:
        dbsession.close()


if __name__ == '__main__':
    main()
