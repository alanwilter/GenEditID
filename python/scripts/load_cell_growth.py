import csv
import sqlalchemy
import sys

from datetime import datetime
from dnascissors.config import cfg
from dnascissors.model import *

def loadIncucyte(plateId, fileName):

    plate = session.query(Plate).filter(Plate.geid == plateId).first()

    if not plate:
        raise Exception('No plate %s' % plateId)

    print("Loading Incucyte growth information for plate %s from %s" % (plateId, fileName))

    with open(fileName, 'r') as icwFile:
        reader = csv.reader(icwFile, delimiter = '\t')

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
                        raise Exception("There is no well at %s%d on plate %s" % (wellrow, wellcolumn, plate.geid))

                    value = float(line[column])

                    print("Setting well %s%s on %s for Incucyte %s at %s" % (well.row, well.column, plate.geid, value, line[0]))

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


if len(sys.argv) < 3:
    print("Need the plate id and at least one Incucyte file.", file=sys.stderr)
    sys.exit(1)

plateId = sys.argv[1]

engine = sqlalchemy.create_engine(cfg['DATABASE_URI'])

Base.metadata.bind = engine

DBSession = sqlalchemy.orm.sessionmaker(bind=engine)

session = DBSession()

for fileIndex in range(2, len(sys.argv)):
    try:
        loadIncucyte(plateId, sys.argv[fileIndex])
    except Exception as e:
        print(e, file=sys.stderr)

session.close()
