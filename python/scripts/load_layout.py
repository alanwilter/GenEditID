import os
import argparse
import log as logger

from geneditid.config import cfg
from geneditid.connect import dbsession
from geneditid.loader import ProjectDataLoader


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--layout", dest="file_layout", action="store", help="The file layout e.g. '20170118_GEP00001.xlsx'", required=True)
    parser.add_argument("--geid", dest="geid", action="store", help="Project GEP ID e.g. 'GEP00001'", required=True)
    options = parser.parse_args()

    log = logger.get_custom_logger(os.path.join(cfg['PROJECTS_FOLDER'], 'load_layout.log'))

    try:
        loader = ProjectDataLoader(dbsession, options.geid, options.file_layout)
        loader.load()
        session.commit()
    except Exception as e:
        log.exception(e)
        session.rollback()
    finally:
        session.close()


if __name__ == '__main__':
    main()
