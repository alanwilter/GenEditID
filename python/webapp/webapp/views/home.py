from pyramid.response import Response
from pyramid.view import view_config

from sqlalchemy.exc import DBAPIError

from dnascissors.model import Project


@view_config(route_name='home', renderer='../templates/home.pt')
def my_view(request):
    try:
        headers = ["geid",
                   "name",
                   "scientist",
                   "start date",
                   "end date",
                   "description",
                   "comments",
                   "nb exp. layouts"]
        rows = []
        projects = request.dbsession.query(Project).all()
        for project in projects:
            rows.append([project.geid,
                         project.name,
                         project.scientist,
                         project.start_date,
                         project.end_date,
                         project.description,
                         project.comments,
                         len(project.experiment_layouts)])
        return dict(title='Genome Editing Core',
                    nb_projects=len(projects),
                    column_headers=headers,
                    rows=rows)
    except DBAPIError:
        return Response(db_err_msg, content_type='text/plain', status=500)


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
