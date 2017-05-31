import colander
import deform.widget

from pyramid.httpexceptions import HTTPFound
from pyramid.view import view_config

from dnascissors.model import Project

# See http://docs.pylonsproject.org/projects/pyramid/en/latest/quick_tutorial/forms.html

# File uploads: http://docs.pylonsproject.org/projects/pyramid-cookbook/en/latest/forms/file_uploads.html

'''
class SelectContent(colander.MappingSchema):
    id = colandar.SchemaNode(colander.Integer())
    name = colandar.SchemaNode(colander.String())
'''

class ProjectViews(object):
    
    def __init__(self, request):
        self.request = request
        self.dbsession = request.dbsession
    
    @property
    def projects_form(self):
        schema = SelectContent()
        return deform.Form(schema, buttons=('submit',))

    @view_config(route_name="projects", renderer="../templates/selectproject.pt")
    def view_projects(self):
        return dict(projects=self.dbsession.query(Project).all(), title="Gene Editing Projects")
    
    @view_config(route_name="project_view", renderer="../templates/viewproject.pt")
    def view_project(self):
        id = self.request.matchdict['projectid']
        project = self.dbsession.query(Project).filter(Project.id == id).first()
        return dict(project=project, title="Gene Editing Project %s" % project.geid)
    
'''    
    @view_config(route_name="projects", renderer="../templates/selectproject.pt")
    def view_projects(self):
        
        form = self.projects_form().render()
        
        projects = self.dbsession.query(Project).all()
        map = dict()
        for p in projects:
            map[p.id] = p.name
        return dict(projects_form=form, projects=map)
'''
    