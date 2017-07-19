import os
import uuid
import shutil
import logging

from pyramid.response import Response
from pyramid.view import view_config
from pyramid.httpexceptions import HTTPFound

from sqlalchemy.exc import DBAPIError

from dnascissors.model import Project
from dnascissors.loader import ExistingEntityException
from dnascissors.loader import LayoutLoader


class HomeViews(object):
    def __init__(self, request):
        self.logger = logging.getLogger('webapp')
        self.request = request
        self.dbsession = request.dbsession

    @view_config(route_name='home', renderer='../templates/home.pt')
    def ge_home_page(self):
        return_map = {'title': 'Genome Editing Core',
                      'clash': False,
                      'error': False,
                      'column_headers': ["geid",
                                         "project details",
                                         "name",
                                         "scientist",
                                         "group leader",
                                         "group",
                                         "start date",
                                         "end date",
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
            for project in projects:
                rows.append([project.geid,
                             "project/{}".format(project.id),
                             project.name,
                             project.scientist,
                             project.group_leader,
                             project.group,
                             project.start_date,
                             project.end_date,
                             project.description,
                             project.comments,
                             project.is_abundance_data_available,
                             project.is_growth_data_available,
                             project.is_variant_data_available,
                             ])
            return_map['rows'] = rows
            # Create new project
            if 'submit' in self.request.params:
                clean_existing = False
                try:
                    clean_existing = self.request.POST['blat']
                except KeyError:
                    pass
                file_path = None
                try:
                    file_path = self._upload("layoutfile", ".xlsx")
                    # Now do the load, after cleaning project.
                    loader = LayoutLoader(self.dbsession, file_path)
                    try:
                        loader.load_all(clean_existing)
                        url = self.request.route_url('home')
                        return HTTPFound(url)
                    except ExistingEntityException as e:
                        return_map['error'] = e.message
                        return_map['clash'] = e
                        return_map['layoutfile'] = self.request.POST['layoutfile']
                except Exception as e:
                    self.logger.error("Have an unexpected error while creating project: {}".format(e))
                    return_map['error'] = str(e)
                finally:
                    if file_path:
                        try:
                            os.remove(file_path)
                        except OSError:
                            pass
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

    def _upload(self, property, suffix):
        filename = self.request.POST[property].filename
        filedata = self.request.POST[property].file
        if not filedata:
            return None
        self.logger.debug("Uploaded = %s" % filename)
        file_path = os.path.join('uploads/', "{}{}".format(uuid.uuid4(), suffix))
        temp_file_path = file_path + '~'
        try:
            filedata.seek(0)
            with open(temp_file_path, 'wb') as output_file:
                shutil.copyfileobj(filedata, output_file)
            os.rename(temp_file_path, file_path)
            statinfo = os.stat(file_path)
            self.logger.info("Uploaded a file of {:d} bytes".format(statinfo.st_size))
        finally:
            try:
                os.remove(temp_file_path)
            except OSError:
                pass
        return file_path
