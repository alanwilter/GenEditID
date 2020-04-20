import os
import argparse
import sqlalchemy

import log as logger
from geneditid.config import cfg
from geneditid.model import Base
from geneditid.plotter import Plotter


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--geid", dest="geid", action="store", help="Project GEP ID e.g. 'GEP00001'", required=True)
    options = parser.parse_args()

    log = logger.get_custom_logger(os.path.join(cfg['PROJECTS_FOLDER'], 'run_variandid.log'))

    engine = sqlalchemy.create_engine(cfg['DATABASE_URI'])
    Base.metadata.bind = engine
    DBSession = sqlalchemy.orm.sessionmaker(bind=engine)
    session = DBSession()
    try:
        plotter = Plotter(session, options.geid)
        plotter.coverage_plot()
        plotter.variant_impact_plot()
        #plotter.heatmap_plot()

    except Exception as e:
        logging.exception(e)
    finally:
        session.close()

if __name__ == '__main__':
    main()
