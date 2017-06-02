import colander
import deform.widget

from pyramid.httpexceptions import HTTPFound
from pyramid.view import view_config

from dnascissors.model import Project

# See http://docs.pylonsproject.org/projects/pyramid/en/latest/quick_tutorial/forms.html

# File uploads: http://docs.pylonsproject.org/projects/pyramid-cookbook/en/latest/forms/file_uploads.html

class ProjectContent(colander.MappingSchema):
    # id = colander.SchemaNode(colander.Integer())
    geid = colander.SchemaNode(colander.String(), title="GEID")
    name = colander.SchemaNode(colander.String(), title="Name")

class ProjectViews(object):
    
    def __init__(self, request):
        self.request = request
        self.dbsession = request.dbsession
    
    @property
    def projects_form(self):
        schema = ProjectContent().bind(request=self.request)
        return deform.Form(schema, buttons=('submit',))

    @view_config(route_name="projects", renderer="../templates/selectproject.pt")
    def view_projects(self):
        return dict(projects=self.dbsession.query(Project).all(), title="Gene Editing Projects")
    
    @view_config(route_name="project_view", renderer="../templates/viewproject.pt")
    def view_project(self):
        id = self.request.matchdict['projectid']
        project = self.dbsession.query(Project).filter(Project.id == id).first()
        return dict(project=project, title="Gene Editing Project %s" % project.geid)
    
    @view_config(route_name="project_edit", renderer="../templates/editproject.pt")
    def edit_project(self):
        id = self.request.matchdict['projectid']
        project = self.dbsession.query(Project).filter(Project.id == id).first()
        
        title = "Gene Editing Project %s" % project.geid
        
        edit_form = self.projects_form
        
        if 'submit' in self.request.params:
            
            fields = self.request.POST.items()
            
            try:
                appstruct = edit_form.validate(fields)
            except deform.ValidationFailure as e:
                return dict(project=project, form=e.render(), title=title)
            
            print("New id = %s" % appstruct['geid'])
            print("New name = %s" % appstruct['name'])
            
            project.geid = appstruct['geid']
            project.name = appstruct['name']
            
            url = self.request.route_url('project_view', projectid=project.id)
            return HTTPFound(url)
        
        projectMap = dict(id=project.id, geid=project.geid, name=project.name)
        
        form = edit_form.render(projectMap)
        
        return dict(edit_form=form, projectid=project.id, project=project, title=title)

        
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
    