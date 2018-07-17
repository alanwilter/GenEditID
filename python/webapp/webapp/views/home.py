import os
import uuid
import shutil
import logging
import datetime

from pyramid.response import Response
from pyramid.view import view_config
from pyramid.httpexceptions import HTTPFound

from sqlalchemy.exc import DBAPIError

from dnascissors.model import Project


class HomeViews(object):
    def __init__(self, request):
        self.logger = logging.getLogger(__name__)
        self.request = request
        self.dbsession = request.dbsession

    @view_config(route_name='home', renderer='../templates/home.pt')
    def ge_home_page(self):
        return_map = {'title': 'Genome Editing Core',
                      'clash': False,
                      'error': False,
                      'info': False,
                      'column_headers': ["geid",
                                         "view",
                                         "edit",
                                         "sequence",
                                         "name",
                                         "scientist",
                                         "group leader",
                                         "group",
                                         "date",
                                         "description",
                                         "comments",
                                         "abundance data",
                                         "growth data",
                                         "ngs data",
                                         ]
                        }
        try:
            # Table of projects
            rows = []
            projects = self.dbsession.query(Project).all()
            last_project = self.dbsession.query(Project).order_by(Project.id.desc()).first()
            for project in projects:
                rows.append([project.geid,
                             "project/{}".format(project.id),
                             "project/{}/edit".format(project.id),
                             "project/{}/sequence".format(project.id),
                             project.name,
                             project.scientist,
                             project.group_leader,
                             project.group,
                             project.start_date,
                             project.description,
                             project.comments,
                             project.is_abundance_data_available,
                             project.is_growth_data_available,
                             project.is_variant_data_available,
                             ])
            return_map['rows'] = rows
            # New project creation form
            if 'submit_project' in self.request.params:
                fields = dict(self.request.POST.items())
                try:
                    next_geid = "GEP{:05d}".format(int(last_project.geid[3:]) + 1)
                    self.dbsession.begin_nested()
                    project = Project()
                    project.geid = next_geid
                    project.name = fields['project_name']
                    project.project_type = fields['project_type']
                    project.start_date = datetime.date.today()
                    project.scientist = fields['project_scientist']
                    project.group = fields['project_group']
                    project.group_leader = fields['project_group_leader']
                    project.description = fields['project_description'].strip()[:1024]
                    self.dbsession.add(project)
                    self.dbsession.commit()
                    return_map['info'] = "Project {} has been created.".format(project.geid)
                    return HTTPFound(self.request.route_url('home'))
                except ValueError as e:
                    self.dbsession.rollback()
                    self.logger.error("Unexpected value error: {}".format(e))
                    return_map['error'] = str(e)
                except DBAPIError as e:
                    self.dbsession.rollback()
                    self.logger.error("Unexpected database error: {}".format(e))
                    return_map['error'] = str(e)
                finally:
                    rows = []
                    projects = self.dbsession.query(Project).all()
                    for project in projects:
                        rows.append([project.geid,
                                     "project/{}".format(project.id),
                                     "project/{}/edit".format(project.id),
                                     "project/{}/sequence".format(project.id),
                                     project.name,
                                     project.scientist,
                                     project.group_leader,
                                     project.group,
                                     project.start_date,
                                     project.description,
                                     project.comments,
                                     project.is_abundance_data_available,
                                     project.is_growth_data_available,
                                     project.is_variant_data_available,
                                     ])
                    return_map['rows'] = rows
                    return return_map
            return return_map
        except DBAPIError:
            return Response(self.db_err_msg, content_type='text/plain', status=500)

    db_err_msg = """\
    Pyramid is having a problem using your SQL database.  The problem
    might be caused by one of the following things:

    1.  You may need to run the "initialize_webapp_db" script
        to initialize your database tables.  Check your virtual
        environment's "bin" directory for this script and try to run it.

    2.  Your database server may not be running.  Check that the
        database server referred to by the "sqlalchemy.url" setting in
        your "development.ini" file is running.

    After you fix the problem, please restart the Pyramid application to
    try it again.
    """
