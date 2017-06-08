import colander
import deform.widget
import datetime

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
    
    def projects_form(self, buttonTitle):
        schema = ProjectContent().bind(request=self.request)
        submitButton = deform.form.Button(name='submit', title=buttonTitle)
        return deform.Form(schema, buttons=(submitButton,))

    @view_config(route_name="projects", renderer="../templates/project/selectproject.pt")
    def view_projects(self):
        return dict(projects=self.dbsession.query(Project).all(), title="Gene Editing Projects")
    
    @view_config(route_name="project_view", renderer="../templates/project/viewproject.pt")
    def view_project(self):
        id = self.request.matchdict['projectid']
        project = self.dbsession.query(Project).filter(Project.id == id).one()
        return dict(project=project, title="Gene Editing Project %s" % project.geid)
    
    @view_config(route_name="project_add", renderer="../templates/project/addproject.pt")
    def add_project(self):
        title = "Create New Gene Editing Project"
        
        add_form = self.projects_form("Create Project")
        
        if 'submit' in self.request.params:
            
            fields = self.request.POST.items()
            
            try:
                appstruct = add_form.validate(fields)
            except deform.ValidationFailure as e:
                return dict(project=project, form=e.render(), title=title)
            
            print("New id = %s" % appstruct['geid'])
            print("New name = %s" % appstruct['name'])
            
            project = Project()
            project.geid = appstruct['geid']
            project.name = appstruct['name']
            project.start_date = datetime.datetime.date(datetime.datetime.now())
            
            self.dbsession.add(project)
            
            url = self.request.route_url('projects')
            return HTTPFound(url)
        
        return dict(title=title)

    @view_config(route_name="project_edit", renderer="../templates/project/editproject.pt")
    def edit_project(self):
        id = self.request.matchdict['projectid']
        project = self.dbsession.query(Project).filter(Project.id == id).one()
        
        title = "Gene Editing Project %s" % project.geid
        
        edit_form = self.projects_form("Update")
        
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
        
        return dict(projectid=project.id, project=project, title=title)
