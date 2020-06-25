import logging
import datetime
import os

from pyramid.response import Response
from pyramid.view import view_config
from pyramid.httpexceptions import HTTPFound

from sqlalchemy.exc import DBAPIError

from geneditid.config import cfg
from geneditid.model import Project


class HomeViews(object):
    def __init__(self, request):
        self.logger = logging.getLogger(__name__)
        self.request = request
        self.dbsession = request.dbsession
        self.projects_folder = cfg['PROJECTS_FOLDER']

    @view_config(route_name='home', renderer='../templates/home.pt')
    def ge_home_page(self):
        return_map = {'title': 'GenEditID',
                      'error': False,
                      'info': False,
                      'column_headers': ["project identifier",
                                         "project link",
                                         "name",
                                         "creation date",
                                         "data location",
                                         "sequencing data",
                                         "abundance data",
                                         "growth data",
                                         ]
                        }
        try:
            # Table of projects
            rows = []
            projects = self.dbsession.query(Project).all()
            last_project = self.dbsession.query(Project).order_by(Project.id.desc()).first()
            for project in projects:
                rows.append([project.geid,
                             "project/{}".format(project.geid),
                             project.name,
                             project.start_date,
                             project.project_folder,
                             project.is_sequencing_data_available,
                             project.is_abundance_data_available,
                             project.is_growth_data_available,
                             ])
            return_map['rows'] = rows
            # New project creation form
            if 'submit_project' in self.request.params:
                fields = dict(self.request.POST.items())
                self.logger.debug(fields)
                try:
                    if last_project:
                        next_geid = "GEP{:05d}".format(int(last_project.geid[3:]) + 1)
                    else:
                        next_geid = "GEP00001"
                    self.dbsession.begin_nested()
                    project = Project()
                    project.geid = next_geid
                    project.name = fields['project_name']
                    project.start_date = datetime.date.today()
                    project.scientist = fields['project_scientist']
                    project.group = fields['project_group']
                    project.group_leader = fields['project_group_leader']
                    project.description = fields['project_description'].strip()[:1024]
                    self.dbsession.add(project)
                    self.dbsession.commit()
                    if not os.path.exists(project.project_folder):
                        os.makedirs(project.project_folder)
                        os.makedirs(os.path.join(project.project_folder, cfg['FASTQ_SUBFOLDER']))
                    return_map['info'] = "Project {} has been created.".format(project.geid)
                    return HTTPFound(self.request.route_url('home'))
                except ValueError as e:
                    self.dbsession.rollback()
                    self.logger.error("Unexpected value error: {}".format(e))
                    return_map['error'] = "Unexpected value error: {}".format(e)
                except DBAPIError as e:
                    self.dbsession.rollback()
                    self.logger.error("Unexpected database error: {}".format(e))
                    if 'UNIQUE constraint failed: project.name' in str(e):
                        return_map['error'] = 'Project name "{}" already exists, please choose another name.'.format(fields['project_name'])
                    else:
                        return_map['error'] = "Unexpected database error: {}".format(e)
                finally:
                    rows = []
                    projects = self.dbsession.query(Project).all()
                    for project in projects:
                        rows.append([project.geid,
                                     "project/{}".format(project.geid),
                                     project.name,
                                     project.start_date,
                                     project.project_folder,
                                     project.is_sequencing_data_available,
                                     project.is_abundance_data_available,
                                     project.is_growth_data_available,
                                     ])
                    return_map['rows'] = rows
                    return return_map
            return return_map
        except DBAPIError:
            return Response(self.db_err_msg, content_type='text/plain', status=500)

    db_err_msg = """\
    Pyramid is having a problem using your SQL database.  The problem
    might be caused by one of the following things:

    1.  You may need to run the "initialize_geneditidapp_db" script
        to initialize your database tables.  Check your virtual
        environment's "bin" directory for this script and try to run it.

    2.  Your database server may not be running.  Check that the
        database server referred to by the "sqlalchemy.url" setting in
        your "development.ini" file is running.

    After you fix the problem, please restart the Pyramid application to
    try it again.
    """
