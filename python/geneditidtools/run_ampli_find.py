import os
import argparse

import geneditid.log as logger
from geneditid.config import cfg
from geneditid.connect import dbsession
from geneditid.finder import AmpliconFinder


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--geid", dest="geid", action="store", help="Project GEP ID e.g. 'GEP00001'", required=True)
    options = parser.parse_args()

    log = logger.get_custom_logger(os.path.join(cfg['PROJECTS_FOLDER'], 'run_ampli_find.log'))

    try:
        log.info('Searching amplicon sequences for project ID\t{}'.format(options.geid))
        finder = AmpliconFinder(dbsession, options.geid)
        finder.write_amplicount_config_file()
    except Exception as e:
        logging.exception(e)
    finally:
        dbsession.close()


if __name__ == '__main__':
    main()
