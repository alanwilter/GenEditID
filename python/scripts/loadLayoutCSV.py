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
    selectionFile = None
    layoutFile = None
    
    def toNone(self, str):
        if str == None:
            return None
        
        str = str.strip()
        
        if str == '':
            return None
        
        return str
    
    def toStrand(self, str, lineNumber):
        str = self.toNone(str)
        
        if not str:
            return None
        
        if str in ['+', 'positive', 'p', 'forward', 'f']:
            return 'forward'
        
        if str in ['-', 'negative', 'n', 'reverse', 'r']:
            return 'reverse'
        
        raise ValueError('Strand must be "forward" or "reverse", not "%s" (line %d)' % (str, lineNumber))

    
    def loadAll(self, session):
        self.loadProjects(session)
        session.commit()
        self.loadTargets(session)
        session.commit()
        self.loadGuides(session)
        session.commit()
        self.loadAmpliconSelection(session)
        session.commit()

    
    def loadProjects(self, session):
        with open(self.projectsFile, 'r') as fileReader:
            reader = csv.reader(fileReader, delimiter=self.columnSeparator)

            headers = next(reader)
            
            lineNumber = 1
            
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
            
            lineNumber = 1
            
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
                target.strand = self.toStrand(line[8], lineNumber)
                    
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


    def loadAmpliconSelection(self, session):
        with open(self.selectionFile, 'r') as fileReader:
            reader = csv.reader(fileReader, delimiter=self.columnSeparator)
            
            headers = next(reader)
            
            lineNumber = 1
            previous_strand = None
            
            for line in reader:
                
                ++lineNumber
                
                # Need the guide first.
                
                guideName = self.toNone(line[0])
                
                if not guideName:
                    raise Exception('Guide name is required on line %d' % lineNumber)
                
                guide = session.query(Guide).filter(Guide.name == guideName).first()
                
                if not guide:
                    raise Exception('Guide "%s" does not exist (line %d)' % (projectId, lineNumber))
                
                
                # Find or create the amplicon
                
                ampliconName = self.toNone(line[5])
                
                if not ampliconName:
                    raise Exception('Amplicon name is required on line %d' % lineNumber)

                amplicon = session.query(Amplicon).filter(Amplicon.name == ampliconName).first()
                
                if not amplicon:
                    amplicon = Amplicon()
                    amplicon.name = ampliconName
                    amplicon.is_on_target = bool(line[6])
                    amplicon.dna_feature = self.toNone(line[7])
                    amplicon.chromosome = self.toNone(line[8])
                    
                    session.add(amplicon);

                    print('Created amplicon %s' % amplicon.name)
                
                
                # Find or create the amplicon selection
                
                guideLocation = int(line[2])
                
                strand = self.toStrand(line[3], lineNumber)
                if strand:
                    previous_strand = strand
                else:
                    if previous_strand:
                        strand = previous_strand
                    else:
                        raise Exception("Have no guide strand on line %d" % lineNumber)
                
                
                selection = session.query(AmpliconSelection).filter(AmpliconSelection.amplicon == amplicon)\
                                                            .filter(AmpliconSelection.guide_location == guideLocation)\
                                                            .filter(AmpliconSelection.guide_strand == strand).first()
                
                if not selection:
                    
                    selection = AmpliconSelection(guide=guide, amplicon=amplicon)
                    selection.experiment_type = self.toNone(line[1])
                    selection.guide_location = guideLocation
                    selection.guide_strand = strand
                    
                    scoreStr = self.toNone(line[4])
                    if scoreStr != None and scoreStr != 'NA':
                        selection.score = int(scoreStr)
                    
                    session.add(selection)
                    
                    print('Created amplicon selection of amplicon %s at %d' % (amplicon.name, selection.guide_location))

                geid = self.toNone(line[9])
                
                if not geid:
                    raise Exception('Primer GEID is required on line %d' % lineNumber)
                
                primer = session.query(Primer).filter(Primer.geid == geid).first()
                
                if not primer:
                    primer = Primer()
                    primer.geid = self.toNone(line[9])
                    primer.sequence = self.toNone(line[10])
                    primer.strand = self.toStrand(line[11], lineNumber)
                    primer.start = int(line[12])
                    primer.end = int(line[13])
                    primer.description = self.toNone(line[14])
                    
                    session.add(primer)
                
                    print('Created primer %s' % primer.geid)
                    
                
                primer.amplicons.append(amplicon)
                
                print('Linked amplicon %s with primer %s' % (amplicon.name, primer.geid))


engine = sqlalchemy.create_engine(cfg['DATABASE_URI'])

Base.metadata.bind = engine

DBSession = sqlalchemy.orm.sessionmaker(bind=engine)

session = DBSession()

deleteCount = session.query(Primer).delete()
print('Deleted %d primers' % deleteCount)

deleteCount = session.query(Amplicon).delete()
print('Deleted %d amplicons' % deleteCount)

deleteCount = session.query(Project).delete()
print('Deleted %d projects' % deleteCount)

# This is a temporary hack.

datadir = '../../data/20170127_GEP00001/PlateLayoutCSV/'

loader = LayoutLoader()
loader.projectsFile = datadir + 'project.csv'
loader.targetsFile = datadir + 'targets.csv'
loader.guidesFile = datadir + 'guides.csv'
loader.selectionFile = datadir + 'ampliconselection.csv'
loader.layoutFile = datadir + 'platelayout.csv'

loader.loadAll(session)

session.commit()
