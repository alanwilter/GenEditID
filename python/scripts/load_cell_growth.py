import csv
import sqlalchemy

# logging configuration
import log as logger

from datetime import datetime
from dnascissors.config import cfg
from dnascissors.model import Base, Plate, Well, CellGrowth


def loadIncucyte(log, session, plateId, fileName):

    plate = session.query(Plate).filter(Plate.geid == plateId).first()

    if not plate:
        raise Exception('No plate {:s}'.format(plateId))

    log.info("Loading Incucyte growth information for plate {:s} from {:s}".format(plateId, fileName))

    with open(fileName, 'r') as f:
        reader = csv.reader(f, delimiter='\t')

        start = False
        header = None

        for line in reader:

            if len(line) == 0:
                continue

            if line[0] == 'Date Time':
                print("START")
                start = True
                header = line

            elif start:

                time = datetime.strptime(line[0], '%d/%m/%Y %H:%M:%S')
                hour = int(line[1])

                for column in range(2, len(line)):

                    location = header[column]

                    wellrow = location[0]
                    wellcolumn = int(location[1:])

                    well = session.query(Well)\
                                  .filter(Well.experiment_layout == plate.experiment_layout)\
                                  .filter(Well.row == wellrow)\
                                  .filter(Well.column == wellcolumn)\
                                  .first()

                    if not well:
                        raise Exception("There is no well at {:s}{:d} on plate {:s}".format(wellrow, wellcolumn, plate.geid))

                    value = float(line[column])

                    log.info("Setting well {:s}{:d} on {:s} for Incucyte {:f} at {:s}".format(well.row, well.column, plate.geid, value, line[0]))

                    growth = session.query(CellGrowth)\
                                    .filter(CellGrowth.plate == plate)\
                                    .filter(CellGrowth.well == well)\
                                    .filter(CellGrowth.hours == hour)\
                                    .first()

                    if not growth:
                        growth = CellGrowth(well=well, plate=plate)
                        growth.hours = hour
                        growth.timestamp = time
                        session.add(growth)

                    growth.confluence_percentage = value

    session.commit()


def main():
    import argparse
    import os

    parser = argparse.ArgumentParser()
    parser.add_argument("--plateid", dest="plateid", action="store", help="The plate ID e.g. 'GEP00001_01'", required=True)
    parser.add_argument("--file", dest="file", action="store", help="The Incucyte file e.g. 'GEP00001_01 data/20170127_GEP00001/GEP00001_01_incu.txt'.", required=True)
    options = parser.parse_args()

    log = logger.get_custom_logger(os.path.join(os.path.dirname(__file__), 'load_cell_growth.log'))

    engine = sqlalchemy.create_engine(cfg['DATABASE_URI'])
    Base.metadata.bind = engine
    DBSession = sqlalchemy.orm.sessionmaker(bind=engine)
    session = DBSession()

    try:
        loadIncucyte(log, session, options.plateid, options.file)
    except Exception as e:
        log.exception(e)

    session.close()

if __name__ == '__main__':
    main()
