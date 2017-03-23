import csv
import sqlalchemy

from dnascissors.config import cfg
from dnascissors.model import *

engine = sqlalchemy.create_engine(cfg['DATABASE_URI'])

Base.metadata.bind = engine

DBSession = sqlalchemy.orm.sessionmaker(bind=engine)

session = DBSession()

def loadICW(plateName, fileName):

    plate = session.query(Plate).filter(Plate.name == plateName).first()

    if not plate:
        raise Exception('No plate %s' % plateName)

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
                    well = session.query(Well).filter(Well.plate == plate).filter(Well.row == rowchar).filter(Well.column == column).first()

                    print("Setting well %s%s on %s for ICW %s" % (well.row, well.column, well.plate.name, channel))

                    if not well:
                        raise Exception('No well %s%s in plate' % (rowchar, column, plateName))

                    if channel == 700:
                        well.icw_700 = float(line[column - 1])
                    elif channel == 800:
                        well.icw_800 = float(line[column - 1])

                row = row + 1


loadICW("Plate1", "data/240117_ICW_Plate1.csv")
loadICW("Plate2", "data/240117_ICW_Plate2.csv")
loadICW("Plate3", "data/240117_ICW_Plate3.csv")
loadICW("Plate4", "data/240117_ICW_Plate4.csv")
loadICW("Plate5", "data/240117_ICW_Plate5.csv")
loadICW("Plate6", "data/240117_ICW_Plate6.csv")

session.commit()
