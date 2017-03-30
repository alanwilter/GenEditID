import csv
import sqlalchemy

from dnascissors.config import cfg
from dnascissors.model import *

from datetime import datetime

class LayoutLoader:
    
    columnSeparator = ","
    
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
        loadProjects(session)
        loadTargets(session)
        loadGuides(session)

    
    def loadProjects(self, session):
        with open(projectsFile, 'r') as fileReader:
            reader = csv.reader(fileReader, columnSeparator)

            headers = reader.next()
            
            lineNumber = 0
            
            for line in reader:
                ++lineNumber
                
                projectId = toNone(line[0])
                
                if not projectId:
                    raise Exception('Project identifier is required on line %s', lineNumber)
                
                project = Project()
                project.geid = projectId
                project.name = toNone(line[1])
                project.scientist = toNone(line[2])
                project.group = toNone(line[4])
                project.group_leader = toNone(line[5])
                if not toNone(line[6]):
                    project.start_date = datetime.strptime(line[6], '%%%m%d')
                project.description = toNone(line[7])
                
                session.add(project)

                
    def loadTargets(self, session):
        with open(targetsFile, 'r') as fileReader:
            reader = csv.reader(fileReader, columnSeparator)
            
            headers = reader.next()
            
            lineNumber = 0
            
            for line in reader:
                
                ++lineNumber
                
                projectId = toNone(line[0])
                
                if not projectId:
                    raise Exception('Project identifier is required on line %s', lineNumber)
                
                project = session.query(Project).filter(Project.geid == projectId).first()
                
                if not project:
                    raise Exception('Project "%s" does not exist (line %s', (projectId, lineNumber))
                
                targetName = toNone(line[1])
                
                if not targetName:
                    raise Exception('Target name is required on line %s', lineNumber)
                
                target = Target(project=project)
                target.name = targetName
                target.species = toNone(line[2])
                target.assembly = toNone(line[3])
                target.gene_id = toNone(line[4])
                target.chromosome = toNone(line[5])
                target.start = int(line[6])
                target.end = int(line[7])
                
                if line[8] in ['+', 'positive']:
                    target.strand = '+'
                elif line[8] in ['-', 'negative']:
                    target.strand = '+'
                else:
                    target.strand = '*'
                    
                target.description = line[9]
                
                session.add(target)

    
    def loadGuides(self, session):
        with open(guidesFile, 'r') as fileReader:
            reader = csv.reader(fileReader, columnSeparator)
            
            headers = reader.next()
            
            lineNumber = 0
            
            for line in reader:
                
                ++lineNumber
                
                targetName = toNone(line[0])
                
                if not targetName:
                    raise Exception('Target name is required on line %s', lineNumber)
                
                target = session.query(Target).filter(Target.name == targetName).first()
                
                if not target:
                    raise Exception('Target "%s" does not exist (line %s', (projectId, lineNumber))

                guideName = toNone(line[1])
                
                if not guideName:
                    raise Exception('Guide name is required on line %s', lineNumber)
                
                guide = Guide(target=target)
                guide.name = guideName
                guide.guide_sequence = toNone(line[2])
                guide.pam_sequence = toNone(line[3])
                guide.activity = int(line[4])
                guide.exon = int(line[5])
                guide.nuclease = int(line[6])
                
                session.add(guide)


engine = sqlalchemy.create_engine(cfg['DATABASE_URI'])

Base.metadata.bind = engine

DBSession = sqlalchemy.orm.sessionmaker(bind=engine)

session = DBSession()

deleteCount = session.query(Project).delete()
print('Deleted %s projects' % deleteCount)


loader = LayoutLoader()
loader.projectsFile = 'data/20170127_GEP00001/PlateLayoutCSV/project.csv'
loader.targetsFile = 'data/20170127_GEP00001/PlateLayoutCSV/targets.csv'
loader.guidesFile = 'data/20170127_GEP00001/PlateLayoutCSV/guides.csv'

loader.loadAll()

session.commit()
