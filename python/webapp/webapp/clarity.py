import logging
import datetime

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

from dnascissors.model import Base, Plate
from dnascissors.config import cfg


LOAD_PROJECTS_QUERY = '''
select p.projectid, p.name, p.luid, r.researcherid, r.labid
from project p
inner join researcher r on r.researcherid=p.researcherid
inner join principals on r.researcherid=principals.researcherid
where p.closedate is null
and principals.accountlocked = 'f'
and r.labid <> 1
'''

LOAD_RESEARCHERS_QUERY = '''
select distinct r.researcherid, r.firstname, r.lastname, r.labid
from researcher r
inner join principals p on r.researcherid=p.researcherid
where p.accountlocked = 'f'
and r.labid <> 1
'''

LOAD_LABS_QUERY = '''
select l.labid, l.name as lname from lab l
where l.labid <> 1
'''

LOAD_PROJECTS_FOR_USER_QUERY = '''
select p.projectid, p.name as pname, p.luid
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

    def does_slx_exist(self, slx):
        return len(self.api.list_filter('sample', 'udf.SLX Identifier', slx).sample) > 0

    def does_project_exist(self, name):
        return len(self.api.list_filter_by_name('project', name).project) > 0

    def get_lab_researcher_project_map(self):
        map = dict()
        self.db.execute(LOAD_LABS_QUERY)
        for results in self.db.fetchall():
            labid = results['labid']
            values = {'id': labid, 'name': results['lname'], 'researchers': dict()}
            map[labid] = values
        self.db.execute(LOAD_RESEARCHERS_QUERY)
        for results in self.db.fetchall():
            labid = results['labid']
            researcherid = results['researcherid']
            values = {'id': researcherid, 'projects': dict(),
                      'name': "{} {}".format(results['firstname'], results['lastname'])}
            lab = map[labid]
            if lab:
                lab['researchers'][researcherid] = values
        self.db.execute(LOAD_PROJECTS_QUERY)
        for results in self.db.fetchall():
            labid = results['labid']
            researcherid = results['researcherid']
            projectid = results['luid']
            values = {'id': projectid, 'name': results['name']}
            lab = map[labid]
            if lab:
                researcher = lab['researchers'][researcherid]
                if researcher:
                    researcher['projects'][projectid] = values
        return map

    def create_project(self, researcher_id, project_name):
        new_project = glsapi.project.project()
        new_project.name = project_name
        new_project.open_date = "{:%Y-%m-%d}".format(datetime.date.today())
        new_project.researcher = glsapi.project.researcher(uri=self.api.load('researcher', researcher_id).uri)
        # new_project.field.append(glsapi.userdefined.field('5352', name = 'Redmine Issue'))
        new_project = self.api.create('project', new_project)
        return new_project

    def submit_samples(self, project_id, plates, override_slx_id=None):

        # This is faulty for more than one plate. Need to discuss how
        # more than one plate is submitted.
        pool_wells_by_plate = dict()
        slx_id = None
        for plate in plates:
            pool_wells = []
            pool_wells_by_plate[plate.id] = pool_wells
            for well in plate.experiment_layout.wells:
                content = well.well_content
                if content and (content.content_type != 'empty'):
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
            self.log.warn("Overriding database SLX id with {}", override_slx_id)

        existing_pools = self.api.list_filter_by_name('artifact', slx_id).artifact
        if existing_pools:
            ep = existing_pools[0]
            raise Exception("Already have a pool for {}: {}".format(slx_id, ep.limsid))

        project = self.api.load('project', project_id)

        container_type = self.api.list_filter_by_name('container_type', '96 well plate').container_type[0]

        containers_by_plate = dict()
        for plate in plates:
            # Reuse an empty container with the same name if one exists.
            container_search = { 'name':plate.barcode, 'state':'Empty', 'type':container_type.name }
    
            empty_containers = self.api.list_filters('container', container_search)

            if empty_containers.container:
                container = self.api.load('container', empty_containers.container[0].limsid)
                self.log.info("Reusing empty {} {} for {}".format(container_type.name, container.limsid, plate.barcode))
                containers_by_plate[plate.id] = container
            else:
                container = glsapi.container.container()
                container.type = glsapi.container.container_type(uri=container_type.uri)
                container.name = plate.barcode
                container = self.api.create('container', container)
                self.log.info("Created new {} {} for {}".format(container_type.name, container.limsid, plate.barcode))
                containers_by_plate[plate.id] = container

        samples = []
        artifacts = []

        try:
            for plate in plates:
                pool_wells = pool_wells_by_plate[plate.id]
                container = containers_by_plate[plate.id]
                for well in pool_wells:
                    for sl_content in well.sequencing_library_contents:
                        sl_lib = sl_content.sequencing_library
                        sample = glsapi.sample.samplecreation()
                        sample.location = glsapi.ri.location(container = glsapi.ri.container(uri=container.uri), value_="{}:{}".format(well.row, well.column))
                        sample.project = glsapi.sample.project(uri=project.uri)
                        sample.name = sl_content.sequencing_sample_name
                        sample.field.append(glsapi.userdefined.field('DNA', name='Sample Type'))
                        sample.field.append(glsapi.userdefined.field('Cell Line', name='Sample Source'))
                        sample.field.append(glsapi.userdefined.field('Amplicon Low-Diversity', name='Library Type'))
                        sample.field.append(glsapi.userdefined.field('Homo sapiens', name='Reference Genome'))
                        sample.field.append(glsapi.userdefined.field(slx_id, name='SLX Identifier'))
                        sample.field.append(glsapi.userdefined.field('MiSeq Nano', name='Workflow'))
                        sample.field.append(glsapi.userdefined.field('Paired End', name='Sequencing Type'))
                        sample.field.append(glsapi.userdefined.field('250', name='Read Length'))
                        sample.field.append(glsapi.userdefined.field(sl_lib.library_type, name='Index Type'))
                        sample.field.append(glsapi.userdefined.field('1', name='Number of Lanes'))
                        sample.field.append(glsapi.userdefined.field('Cell Line', name='Sample Source'))
                        sample.field.append(glsapi.userdefined.field('289', name='Average Library Length'))  # needs to go to UI
                        sample.field.append(glsapi.userdefined.field('SWAG/076', name='Billing Information'))  # needs to get it from project
                        sample.field.append(glsapi.userdefined.field('Standard', name='Priority Status'))
                        sample.field.append(glsapi.userdefined.field('Please use 20% PhiX.', name='Submission Comments'))
                        sample.field.append(glsapi.userdefined.field('SLX Version 44', name='Version Number'))
                        sample.field.append(glsapi.userdefined.field(well.row, name='Row'))
                        sample.field.append(glsapi.userdefined.field(well.column, name='Column'))
                        sample.field.append(glsapi.userdefined.field('-1', name='Concentration'))
                        sample.field.append(glsapi.userdefined.field('-1', name='Volume'))
                        sample.field.append(glsapi.userdefined.field('Not Assigned', name='Sequencer'))
                        sample.field.append(glsapi.userdefined.field(len(pool_wells), name='Pool Size'))
                        sample = self.api.create('sample', sample)
                        samples.append(sample)
                        self.log.debug("New sample id = {}".format(sample.limsid))
    
                        artifact = self.api.load_by_uri('artifact', sample.artifact.uri)
                        artifact.reagent_label.append(glsapi.artifact.reagent_label(name=sl_content.sequencing_barcode))
                        artifact = self.api.update(artifact)
                        artifacts.append(artifact)
                        self.log.debug("Sample {} barcode set to {}".format(sample.limsid, artifact.reagent_label[0].name))

        except Exception as e:
            self.log.error(e)
            self.log.error("Creating samples has failed. Need to remove existing containers.")
            for container in containers_by_plate.values():
                self.api.delete(container)
            raise e

        # If all have been created, route them into MiSeq Express Nano work flow
        # <workflow status="ACTIVE" uri="https://limsdev.cruk.cam.ac.uk/api/v2/configuration/workflows/1452" name="SLX: Sequencing for MiSeq Express Nano v1"/>
        # Anne: I would not route the samples, but let the automatic assigment script does it
        # routing = glsapi.routing.routing()
        # assignment = glsapi.routing.extArtifactAssignments(workflow_uri="https://limsdev.cruk.cam.ac.uk/api/v2/configuration/workflows/1452")
        # routing.assign.append(assignment)
        # for a in artifacts:
        #     assignment.artifact.append(glsapi.routing.artifact(uri=a.uri))
        #
        # self.api.create('routing', routing)
        # self.log.info("Samples routed into {}".format(routing.assign[0].workflow_uri))


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
    claire_id = 80

    # Get projects for me (Rich).
    for infomap in clarity.get_projects_for_researcher(rich_id).values():
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
