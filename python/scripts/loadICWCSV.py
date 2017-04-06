import csv
import sqlalchemy
import sys

from dnascissors.config import cfg
from dnascissors.model import *

def loadICW(plateId, fileName):

    plate = session.query(Plate).filter(Plate.geid == plateId).first()

    if not plate:
        raise Exception('No plate %s' % plateId)

    print("Loading InCell Western signal for plate %s from %s" % (plateId, fileName))

    with open(fileName, 'r') as icwFile:
        reader = csv.reader(icwFile, delimiter = '\t')

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
                        raise Exception("There is no well at %s%d on plate %s" % (rowchar, column, plate.geid))

                    print("Setting well %s%s on %s for ICW %s" % (well.row, well.column, plate.geid, channel))

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



if len(sys.argv) < 3:
    print("Need the plate id and at least one InCell Western file.", file=sys.stderr)
    sys.exit(1)

plateId = sys.argv[1]

engine = sqlalchemy.create_engine(cfg['DATABASE_URI'])

Base.metadata.bind = engine

DBSession = sqlalchemy.orm.sessionmaker(bind=engine)

session = DBSession()

for fileIndex in range(2, len(sys.argv)):
    try:
        loadICW(plateId, sys.argv[fileIndex])
    except Exception as e:
        print(e, file=sys.stderr)
