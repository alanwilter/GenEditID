import os
import argparse

import geneditid.log as logger
from geneditid.config import cfg
from geneditid.connect import dbsession

from geneditid.plotter import Plotter


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--geid", dest="geid", action="store", help="Project GEP ID e.g. 'GEP00001'", required=True)
    options = parser.parse_args()

    log = logger.get_custom_logger(os.path.join(cfg['PROJECTS_FOLDER'], 'run_variandid.log'))

    try:
        plotter = Plotter(dbsession, options.geid)
        plotter.coverage_plot()
        plotter.variant_impact_plot()
        #plotter.heatmap_plot()
    except Exception as e:
        logging.exception(e)
    finally:
        dbsession.close()

if __name__ == '__main__':
    main()
