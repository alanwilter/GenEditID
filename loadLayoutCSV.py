import csv
import sqlalchemy

from config import cfg
from model import *

engine = sqlalchemy.create_engine(cfg['DATABASE_URI'])

Base.metadata.bind = engine

DBSession = sqlalchemy.orm.sessionmaker(bind=engine)

session = DBSession()

deleteCount = session.query(Project).delete()
print('Deleted %s projects' % deleteCount)

deleteCount = session.query(CellLine).delete()
print('Deleted %s cell lines' % deleteCount)

project = Project(name = 'initial')
session.add(project)

target = Target(name = 'target', project = project)
session.add(target)

oligo = Oligo(name = 'oligo', target = target)
session.add(target)

cellline = CellLine(name = 'cellline')
session.add(cellline)

clonesByName = dict()

def addWell(plate, row, column):
    
    cell = line[column]
    clone = None
    
    if cell != '' and cell != '-':
        if cell in clonesByName:
            clone = clonesByName[cell]
        else:
            clone = Clone(name = cell, oligo = oligo, cell_line = cellline)
            session.add(clone)
            clonesByName[clone.name] = clone

    well = Well(plate = plate, clone = clone, location = '%s%s' % (row, column), empty = (clone == None))
    session.add(well)


with open('data/PlatesLayout_270117.csv', 'r') as layoutFile:
    reader = csv.reader(layoutFile, delimiter = '\t')
    previousLine = []
    
    plate = None
    
    for line in reader:
        
        if line[0] == 'A':
            # This cell indicates the start of a new layout.
            
            plate = Plate(name = previousLine[14])
            session.add(plate)
            
            row = line[0]
            for column in range(1, 12):
                addWell(plate, row, column)
            
        elif line[0] == '':
            # Indicates the end of a table. May also get it at the start.
            
            plate = None
            
        elif plate:
            
            row = line[0]
            for column in range(1, 12):
                addWell(plate, row, column)
        
        previousLine = line
        
session.commit()

