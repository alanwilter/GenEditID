import colander
import deform.widget

from pyramid.httpexceptions import HTTPFound
from pyramid.view import view_config

from dnascissors.model import Project

from webapp.plots.plotter import Plotter


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

    @view_config(route_name="project_view", renderer="../templates/project/viewproject.pt")
    def view_project(self):
        id = self.request.matchdict['projectid']
        project = self.dbsession.query(Project).filter(Project.id == id).one()
        plotter = Plotter()
        layouts = project.experiment_layouts
        headers = ["Plate", "Well", "Sample", "Barcode", "Protein", "Type", "Allele",
                   "Allele Fraction", "Frame", "Variant Type" ]
        rows = []

        for layout in layouts:
            for well in layout.wells:
                for slc in well.sequencing_library_contents:
                    for vc in slc.variant_results:
                        row = []
                        row.append(layout.geid)
                        row.append("{:s}{:02}".format(well.row, well.column))
                        row.append(slc.sequencing_sample_name)
                        row.append(slc.sequencing_barcode)
                        row.append(vc.protein_effect)
                        row.append(vc.consequence)
                        row.append(vc.alleles)
                        row.append("{0:.3f}".format(vc.allele_fraction))
                        row.append(vc.frame)
                        row.append(vc.variant_type)
                        rows.append(row)

        return dict(project=project,
                    title="Gene Editing Project %s" % project.geid,
                    cellgrowthplot=plotter.growth_plot(self.dbsession, project.geid),
                    proteinabundanceplot=plotter.abundance_plot(self.dbsession, project.geid),
                    column_headers=headers,
                    rows=rows)

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
            print("New name = %s" % appstruct['name'])
            project.name = appstruct['name']
            url = self.request.route_url('project_view', projectid=project.id)
            return HTTPFound(url)
        return dict(projectid=project.id, project=project, title=title)
