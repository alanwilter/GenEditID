import sys
import os
import logging

import psycopg2
from psycopg2.extras import RealDictCursor
import sqlalchemy

from glsclient.glsclient import GlsClientApi
from glsclient.config import SERVER, TEST_SERVER, USERNAME, PASSWORD, DB_NAME, FILES_DB_NAME, DB_USERNAME, DB_PASSWORD

import glsapi.artifact
import glsapi.container
import glsapi.project
import glsapi.ri
import glsapi.routing
import glsapi.sample
import glsapi.userdefined

from dnascissors.model import Base, Plate, ExperimentLayout, Well, WellContent
from dnascissors.config import cfg


LOAD_PROJECTS_QUERY = '''
select p.projectid, p.name as pname, r.firstname, r.lastname, l.name as lname
from project p
inner join researcher r on p.researcherid=r.researcherid
inner join lab l on r.labid=l.labid
'''

LOAD_USERS_QUERY = '''
select r.researcherid, r.firstname, r.lastname, l.name as lname
from researcher r
inner join lab l on r.labid=l.labid
'''

LOAD_PROJECTS_FOR_USER_QUERY = '''
select p.projectid, p.name as pname
from project p
where p.researcherid = {}
'''

class ClaritySubmitter(object):
    
    # server=SERVER, username=USERNAME, password=PASSWORD, db_name=DB_NAME, db_username=DB_USERNAME, db_password=DB_PASSWORD
    def __init__(self, use_dev_lims=True):
        self.log = logging.getLogger(__name__)
        lims_server = SERVER
        if use_dev_lims:
            lims_server = TEST_SERVER
        self.api = GlsClientApi(lims_server, USERNAME, PASSWORD)
        self.db_connection = psycopg2.connect(database=DB_NAME, user=DB_USERNAME, password=DB_PASSWORD, host=lims_server, cursor_factory=RealDictCursor)
        self.db = self.db_connection.cursor()
        self.files_db_connection = psycopg2.connect(database=FILES_DB_NAME, user=DB_USERNAME, password=DB_PASSWORD, host=lims_server, cursor_factory=RealDictCursor)
        self.files_db = self.files_db_connection.cursor()

    def close_db_connection(self):
        self.db.close()
        self.db_connection.close()
        self.files_db.close()
        self.files_db_connection.close()

    def get_projects(self):
        result_map = dict()
        self.db.execute(LOAD_PROJECTS_QUERY)
        for results in self.db.fetchall():
            id = results['projectid']
            values = { 'id':id, 'name':results['pname'],
                       'researcher':"{} {}".format(results['firstname'], results['lastname']),
                       'lab':results['lname'] }
            result_map[id] = values
        return result_map

    def get_researchers(self):
        result_map = dict()
        self.db.execute(LOAD_USERS_QUERY)
        for results in self.db.fetchall():
            id = results['researcherid']
            values = { 'id':id,
                       'researcher':"{} {}".format(results['firstname'], results['lastname']),
                       'lab':results['lname'] }
            result_map[id] = values
        return result_map

    def get_projects_for_researchers(self, researcher_id):
        result_map = dict()
        self.db.execute(LOAD_PROJECTS_FOR_USER_QUERY.format(researcher_id))
        for results in self.db.fetchall():
            id = results['projectid']
            values = { 'id':id, 'name':results['pname'] }
            result_map[id] = values
        return result_map
    
    def create_project(self, researcher_id, project_name):
        new_project = glsapi.project.project()
        new_project.name = project_name
        new_project.open_date = '2018-05-09'
        new_project.researcher = glsapi.project.researcher(uri = self.api.load('researcher', researcher_id).uri)
        # new_project.field.append(glsapi.userdefined.field('5352', name = 'Redmine Issue'))
        new_project = self.api.create('project', new_project)
        return new_project
    
    def submit_samples(self, project_id, plate, override_slx_id = None):
        
        pool_wells = []
        slx_id = None
        
        for well in plate.experiment_layout.wells:
            content = well.well_content
            
            if content and content.content_type != 'empty':
                for sl_content in well.sequencing_library_contents:
                    sl_lib = sl_content.sequencing_library
                    if slx_id:
                        if slx_id != sl_lib.slxid:
                            raise Exception("Have a mix of SLX ids in the submission. There must be one only.")
                    else:
                        slx_id = sl_lib.slxid
                    pool_wells.append(well)
        
        if not slx_id:
            raise Exception("Have no SLX id available.")
        
        if override_slx_id:
            slx_id = override_slx_id
            print("Overriding database SLX id with {}", override_slx_id)

        existing_pools = self.api.list_filter_by_name('artifact', slx_id).artifact
        if existing_pools:
            ep = existing_pools[0]
            raise Exception("Already have a pool for {}: {}".format(slx_id, ep.limsid))
        
        project = self.api.load('project', project_id)
        
        container_type = self.api.list_filter_by_name('container_type', '96 well plate').container_type[0]
        
        container = glsapi.container.container()
        container.type = glsapi.container.container_type(uri = container_type.uri)
        container.name = plate.barcode
        container.type_ = glsapi.userdefined.type(name = 'SLX Container')
        container = self.api.create('container', container)

        samples = []
        artifacts = []
        
        try:
            for well in pool_wells:
                for sl_content in well.sequencing_library_contents:
                    sl_lib = sl_content.sequencing_library
                
                    sample = glsapi.sample.samplecreation()
                    sample.location = glsapi.ri.location(container = glsapi.ri.container(uri = container.uri), value_ = "{}:{}".format(well.row, well.column))
                    sample.project = glsapi.sample.project(uri = project.uri)
                    sample.name = sl_content.sequencing_sample_name
                    sample.field.append(glsapi.userdefined.field('DNA', name = 'Sample Type'))
                    sample.field.append(glsapi.userdefined.field('Cell Line', name = 'Sample Source'))
                    sample.field.append(glsapi.userdefined.field('Amplicon Low-Diversity', name = 'Library Type'))
                    sample.field.append(glsapi.userdefined.field('Homo sapiens [GRCh38_hs38d1]', name = 'Reference Genome'))
                    sample.field.append(glsapi.userdefined.field(slx_id, name = 'SLX Identifier'))
                    sample.field.append(glsapi.userdefined.field('MiSeq', name = 'Workflow'))
                    sample.field.append(glsapi.userdefined.field('Paired End', name = 'Sequencing Type'))
                    sample.field.append(glsapi.userdefined.field('300', name = 'Read Length'))
                    sample.field.append(glsapi.userdefined.field(sl_lib.library_type, name = 'Index Type'))
                    sample.field.append(glsapi.userdefined.field('1', name = 'Number of Lanes'))
                    sample.field.append(glsapi.userdefined.field('Cell Line', name = 'Sample Source'))
                    sample.field.append(glsapi.userdefined.field('323', name = 'Average Library Length'))
                    sample.field.append(glsapi.userdefined.field('SWAG/076', name = 'Billing Information'))
                    sample.field.append(glsapi.userdefined.field('Standard', name = 'Priority Status'))
                    sample.field.append(glsapi.userdefined.field('Add 20% PhiX please', name = 'Submission Comments'))
                    sample.field.append(glsapi.userdefined.field('SLX Version 44', name = 'Version Number'))
                    sample.field.append(glsapi.userdefined.field(well.row, name = 'Row'))
                    sample.field.append(glsapi.userdefined.field(well.column, name = 'Column'))
                    sample.field.append(glsapi.userdefined.field('-1', name = 'Concentration'))
                    sample.field.append(glsapi.userdefined.field('-1', name = 'Volume'))
                    sample.field.append(glsapi.userdefined.field('Not Assigned', name = 'Sequencer'))
                    sample.field.append(glsapi.userdefined.field(len(pool_wells), name = 'Pool Size'))
                    sample = self.api.create('sample', sample)
                    samples.append(sample)
                    print("New sample id = {}".format(sample.limsid))
                    
                    artifact = self.api.load_by_uri('artifact', sample.artifact.uri)
                    artifact.reagent_label.append(glsapi.artifact.reagent_label(name = sl_content.sequencing_barcode))
                    artifact = self.api.update(artifact)
                    artifacts.append(artifact)
                    print("Sample {} barcode set to {}".format(sample.limsid, artifact.reagent_label[0].name))

        except Exception as e:
            self.log.error("Creating samples has failed. Need to remove existing container {}.".format(container.limsid))
            self.api.delete(container)
            raise e
        
        # If all have been created, route them into MiSeq Express work flow
        
        routing = glsapi.routing.routing()
        assignment = glsapi.routing.extArtifactAssignments(workflow_uri = "https://limsdev.cruk.cam.ac.uk/api/v2/configuration/workflows/11")
        routing.assign.append(assignment)
        for a in artifacts:
            assignment.artifact.append(glsapi.routing.artifact(uri = a.uri))
        
        self.api.create('routing', routing)
        
        print("Samples routed into {}".format(routing.assign[0].workflow_uri))

def tests(session, clarity):
    print("")
    print("Projects")

    # Listing all projects
    for infomap in clarity.get_projects().values():
        print("{} {}: {} in {}".format(infomap['id'], infomap['name'], infomap['researcher'], infomap['lab']))

    print("")
    print("Researchers")

    # Listing all researchers
    for infomap in clarity.get_researchers().values():
        print("{} - {} in {}".format(infomap['id'], infomap['researcher'], infomap['lab']))
        
    print("")
    print("My projects")

    rich_id = 153

    # Get projects for me (Rich).
    for infomap in clarity.get_projects_for_researchers(rich_id).values():
        print("{} {}".format(infomap['id'], infomap['name']))

    #print("")
    #print("Create a project")
    
    #new_project = clarity.create_project(rich_id, "Created through the API")
    #print(new_project)
    
    print("")
    print("Submit samples to test project")
    
    project_id = 'BOW13101'
    plate_id = 'GEP00007_01_NGS'
    slx_id = 'SLX-15104'
    
    plate = session.query(Plate).filter(Plate.geid == plate_id).first()
    #print("{} {} {}".format(plate.id, plate.geid, plate.barcode))
    
    clarity.submit_samples(project_id, plate, 'SLX-20003')

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                        datefmt='%m-%d %H:%M')
    clarity = ClaritySubmitter(True)
    
    engine = sqlalchemy.create_engine(cfg['DATABASE_URI'])
    Base.metadata.bind = engine
    DBSession = sqlalchemy.orm.sessionmaker(bind=engine)
    session = DBSession()
    try:
        tests(session, clarity)
    except Exception as e:
        logging.exception(e)
    finally:
        session.close()

if __name__ == '__main__':
    main()
