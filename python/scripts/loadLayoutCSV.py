import csv
import sqlalchemy

from dnascissors.config import cfg
from dnascissors.model import *

from datetime import datetime

class LayoutLoader:
    
    columnSeparator = ','
    
    projectsFile = None
    targetsFile = None
    guidesFile = None
    layoutFile = None
    
    def toNone(self, str):
        if str == None:
            return None
        
        str = str.strip()
        
        if str == '':
            return None
        
        return str

    
    def loadAll(self, session):
        self.loadProjects(session)
        session.commit()
        self.loadTargets(session)
        session.commit()
        self.loadGuides(session)
        session.commit()

    
    def loadProjects(self, session):
        with open(self.projectsFile, 'r') as fileReader:
            reader = csv.reader(fileReader, delimiter=self.columnSeparator)

            headers = next(reader)
            
            lineNumber = 0
            
            for line in reader:
                ++lineNumber
                
                projectId = self.toNone(line[0])
                
                if not projectId:
                    raise Exception('Project identifier is required on line %d' % lineNumber)
                
                project = Project()
                project.geid = projectId
                project.name = self.toNone(line[1])
                project.scientist = self.toNone(line[2])
                project.group = self.toNone(line[4])
                project.group_leader = self.toNone(line[5])
                if self.toNone(line[6]):
                    project.start_date = datetime.strptime(line[6], '%Y%m%d')
                    if not project.start_date:
                        raise Exception('Start date %s has not populated' % line[6])
                project.description = self.toNone(line[7])
                
                session.add(project)
                
                print('Created project %s' % project.name)

                
    def loadTargets(self, session):
        with open(self.targetsFile, 'r') as fileReader:
            reader = csv.reader(fileReader, delimiter=self.columnSeparator)
            
            headers = next(reader)
            
            lineNumber = 0
            
            for line in reader:
                
                ++lineNumber
                
                projectId = self.toNone(line[0])
                
                if not projectId:
                    raise Exception('Project identifier is required on line %d' % lineNumber)
                
                project = session.query(Project).filter(Project.geid == projectId).first()
                
                if not project:
                    raise Exception('Project "%s" does not exist (line %d' % (projectId, lineNumber))
                
                targetName = self.toNone(line[1])
                
                if not targetName:
                    raise Exception('Target name is required on line %d' % lineNumber)
                
                target = Target(project=project)
                target.name = targetName
                target.species = self.toNone(line[2])
                target.assembly = self.toNone(line[3])
                target.gene_id = self.toNone(line[4])
                target.chromosome = self.toNone(line[5])
                target.start = int(line[6])
                target.end = int(line[7])
                
                if line[8] in ['+', 'positive', 'p', 'forward', 'f']:
                    target.strand = 'negative'
                elif line[8] in ['-', 'negative', 'n', 'reverse', 'r']:
                    target.strand = 'reverse'
                else:
                    raise ValueError('Strand must be "forward" or "reverse", not "%s" (line %d)' % (line[8], lineNumber))
                    
                target.description = line[9]
                
                session.add(target)
                
                print('Created target %s' % target.name)

    
    def loadGuides(self, session):
        with open(self.guidesFile, 'r') as fileReader:
            reader = csv.reader(fileReader, delimiter=self.columnSeparator)
            
            headers = next(reader)
            
            lineNumber = 0
            
            for line in reader:
                
                ++lineNumber
                
                targetName = self.toNone(line[0])
                
                if not targetName:
                    raise Exception('Target name is required on line %d' % lineNumber)
                
                target = session.query(Target).filter(Target.name == targetName).first()
                
                if not target:
                    raise Exception('Target "%s" does not exist (line %d)' % (projectId, lineNumber))

                guideName = self.toNone(line[1])
                
                if not guideName:
                    raise Exception('Guide name is required on line %d' % lineNumber)
                
                guide = Guide(target=target)
                guide.name = guideName
                guide.guide_sequence = self.toNone(line[2])
                guide.pam_sequence = self.toNone(line[3])
                guide.activity = int(line[4])
                guide.exon = int(line[5])
                guide.nuclease = line[6]
                
                session.add(guide)
                
                print('Created guide %s' % guide.name)


engine = sqlalchemy.create_engine(cfg['DATABASE_URI'])

Base.metadata.bind = engine

DBSession = sqlalchemy.orm.sessionmaker(bind=engine)

session = DBSession()

deleteCount = session.query(Project).delete()
print('Deleted %s projects' % deleteCount)

# This is a temporary hack.

datadir = '../../data/20170127_GEP00001/PlateLayoutCSV/'

loader = LayoutLoader()
loader.projectsFile = datadir + 'project.csv'
loader.targetsFile = datadir + 'targets.csv'
loader.guidesFile = datadir + 'guides.csv'

loader.loadAll(session)

session.commit()
