import os
import argparse

import geneditid.log as logger
from geneditid.config import cfg
from geneditid.connect import dbsession
from geneditid.loader import ProjectLoader

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--geid", dest="geid", action="store", help="Project GEP ID e.g. 'GEP00001'", required=True)
    options = parser.parse_args()

    log = logger.get_custom_logger(os.path.join(cfg['PROJECTS_FOLDER'], 'delete_project.log'))

    try:
        project = ProjectLoader(dbsession)
        project.delete_project(options.geid)
        dbsession.commit()
    except Exception as e:
        log.exception(e)
        dbsession.rollback()
    finally:
        dbsession.close()


if __name__ == '__main__':
    main()
