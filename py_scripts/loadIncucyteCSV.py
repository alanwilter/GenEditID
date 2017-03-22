import csv
import sqlalchemy

from datetime import datetime
from dnascissors.config import cfg
from dnascissors.model import *

engine = sqlalchemy.create_engine(cfg['DATABASE_URI'])

Base.metadata.bind = engine

DBSession = sqlalchemy.orm.sessionmaker(bind=engine)

session = DBSession()

session.query(IncucyteGrowth).delete()

def loadIncucyte(plateName, fileName):

    plate = session.query(Plate).filter(Plate.name == plateName).first()

    if not plate:
        raise Exception('No plate %s' % plateName)

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

                    well = session.query(Well).filter(Well.plate == plate).filter(Well.row == wellrow).filter(Well.column == wellcolumn).first()

                    value = float(line[column])

                    print("Setting well %s%s on %s for Incucyte %s at %s" % (well.row, well.column, well.plate.name, value, line[0]))

                    if not well:
                        raise Exception('No well %s%s in plate' % (rowchar, column, plateName))

                    incucyte = IncucyteGrowth(well = well, hours = hour, timestamp = time, phased_object_confluence = value)
                    session.add(incucyte)


loadIncucyte("Plate1", "data/180117_MCF7_luc-straw_clone3_CRISPR_plate1.txt")
loadIncucyte("Plate2", "data/180117_MCF7_luc-straw_clone3_CRISPR_plate2.txt")
loadIncucyte("Plate3", "data/180117_MCF7_luc-straw_clone3_CRISPR_plate3.txt")
loadIncucyte("Plate4", "data/180117_MCF7_luc-straw_clone3_CRISPR_plate4.txt")
loadIncucyte("Plate5", "data/180117_MCF7_luc-straw_clone3_CRISPR_plate5.txt")
loadIncucyte("Plate6", "data/180117_MCF7_luc-straw_clone3_CRISPR_plate6.txt")

session.commit()
