import csv
import sqlalchemy

# logging configuration
import log as logger

from dnascissors.config import cfg
from dnascissors.model import Base, Plate, Well, ProteinAbundance


def loadICW(log, session, plateId, fileName):

    plate = session.query(Plate).filter(Plate.geid == plateId).first()

    if not plate:
        raise Exception('No plate {:s}'.format(plateId))

    log.info("Loading InCell Western signal for plate {:s} from {:s}".format(plateId, fileName))

    with open(fileName, 'r') as icwFile:
        reader = csv.reader(icwFile, delimiter='\t')

        channel = None
        row = None
        a = ord('A')

        for line in reader:

            if line[0] == '':
                channel = None
                row = None
            elif line[0] == 'Total Channel 700':
                channel = 700
                row = 0
            elif line[0] == 'Total Channel 800':
                channel = 800
                row = 0
            elif channel:
                rowchar = chr(a + row)

                for column in range(1, 13):
                    well = session.query(Well)\
                                  .filter(Well.experiment_layout == plate.experiment_layout)\
                                  .filter(Well.row == rowchar)\
                                  .filter(Well.column == column)\
                                  .first()

                    if not well:
                        raise Exception("There is no well at {:s}{:d} on plate {:s}".format(rowchar, column, plate.geid))

                    log.info("Setting well {:s}{:d} on {:s} for ICW {:d}".format(well.row, well.column, plate.geid, channel))

                    abundance = session.query(ProteinAbundance)\
                                       .filter(ProteinAbundance.well == well)\
                                       .filter(ProteinAbundance.plate == plate)\
                                       .first()

                    if not abundance:
                        abundance = ProteinAbundance(well=well, plate=plate)
                        session.add(abundance)

                    if channel == 700:
                        abundance.intensity_channel_700 = float(line[column - 1])
                    elif channel == 800:
                        abundance.intensity_channel_800 = float(line[column - 1])

                row = row + 1

    session.commit()


def main():
    import argparse
    import os

    parser = argparse.ArgumentParser()
    parser.add_argument("--plateid", dest="plateid", action="store", help="The plate ID e.g. 'GEP00001_01'", required=True)
    parser.add_argument("--file", dest="file", action="store", help="The InCell Western file e.g. 'data/20170127_GEP00001/GEP00001_01_ICW.csv'.", required=True)
    options = parser.parse_args()

    log = logger.get_custom_logger(os.path.join(os.path.dirname(__file__), 'load_protein_abundance.log'))

    engine = sqlalchemy.create_engine(cfg['DATABASE_URI'])
    Base.metadata.bind = engine
    DBSession = sqlalchemy.orm.sessionmaker(bind=engine)
    session = DBSession()

    try:
        loadICW(log, session, options.plateid, options.file)
    except Exception as e:
        log.exception(e)

    session.close()


if __name__ == '__main__':
    main()
