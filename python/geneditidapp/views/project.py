import os
import uuid
import shutil
import socket
import logging
import colander
import requests
import deform.widget

from collections import OrderedDict

from urllib.parse import quote

from pyramid.httpexceptions import HTTPFound
from pyramid.view import view_config

from json import JSONEncoder
from sqlalchemy.exc import DBAPIError

from geneditid.finder import FinderException
from geneditid.loader import ExistingEntityException
from geneditid.loader import LoaderException
from geneditid.loader import ProjectDataLoader
from geneditid.loader import CellGrowthLoader
from geneditid.loader import ProteinAbundanceLoader

from geneditid.model import Project
from geneditid.model import Layout
from geneditid.model import LayoutContent
from geneditid.model import Plate

from geneditid.plotter import Plotter

# See http://docs.pylonsproject.org/projects/pyramid/en/latest/quick_tutorial/forms.html
# File uploads: http://docs.pylonsproject.org/projects/pyramid-cookbook/en/latest/forms/file_uploads.html


class ProjectContent(colander.MappingSchema):
    comments = colander.SchemaNode(colander.String(), title="comments", default='', missing='')


class ProjectViews(object):

    def __init__(self, request):
        self.logger = logging.getLogger(__name__)
        self.request = request
        self.dbsession = request.dbsession

    def projects_form(self, buttonTitle):
        schema = ProjectContent().bind(request=self.request)
        submitButton = deform.form.Button(name='submit_comments', title=buttonTitle)
        return deform.Form(schema, buttons=(submitButton,))

    def get_project_table(self, project):
        headers = [
            "project identifier",
            "name",
            "type",
            "creation date",
            "data location",
            "sequencing data",
            "abundance data",
            "growth data",
            "scientist",
            "group",
            "description",
            "comments",
            ]
        rows = [[project.geid,
                 project.name,
                 project.project_type,
                 project.start_date,
                 project.project_folder,
                 project.is_sequencing_data_available,
                 project.is_abundance_data_available,
                 project.is_growth_data_available,
                 project.scientist,
                 project.group,
                 project.description,
                 project.comments
                 ]]
        return headers, rows

    def get_target_table(self, project):
        headers = [
            "name",
            "species",
            "assembly",
            "gene",
            "chromosome",
            "start",
            "end",
            "strand",
            "description"]
        rows = []
        for target in project.targets:
            row = []
            row.append(target.name)
            row.append(target.genome.species)
            row.append(target.genome.assembly)
            row.append(target.gene_id)
            row.append(target.chromosome)
            row.append(target.start)
            row.append(target.end)
            row.append(target.strand)
            row.append(target.description)
            rows.append(row)
        return headers, rows

    def get_guide_table(self, project):
        headers = [
            "target name",
            "guide name",
            "guide sequence",
            "pam sequence",
            "activity",
            "exon",
            "nuclease"
        ]
        rows = []
        for guide in project.guides:
            row = []
            row.append(guide.target.name)
            row.append(guide.name)
            row.append(guide.guide_sequence)
            row.append(guide.pam_sequence)
            row.append(guide.activity)
            row.append(guide.exon)
            row.append(guide.nuclease)
            rows.append(row)
        return headers, rows

    def get_amplicon_table(self, project):
        headers = [
            "amplicon name",
            "forward primer",
            "reverse primer",
            "amplicon coordinates",
            "guide name",
            "experiment type",
            "guide location",
            "is on target?",
            "DNA feature",
            "chromosome",
            "score",
            "description"
        ]
        rows = []
        for amplicon in project.amplicons:
            row = []
            row.append(amplicon.name)
            row.append(amplicon.fprimer.sequence)
            row.append(amplicon.rprimer.sequence)
            row.append(amplicon.coordinates)
            row.append(amplicon.guide.name)
            row.append(amplicon.experiment_type)
            row.append(amplicon.guide_location)
            row.append(amplicon.is_on_target)
            row.append(amplicon.dna_feature)
            row.append(amplicon.chromosome)
            row.append(amplicon.score)
            row.append(amplicon.description)
            rows.append(row)
        return headers, rows

    def get_layout_table(self, project):
        headers = [
            "layout id",
            "well position",
            "guide name",
            "sequencing barcode",
            "sequencing dna source",
            "sequencing project id",
            "sequencing sample name",
            "sequencing library type",
            "cell line name",
            "clone name",
            "cell pool",
            "content type",
            "is control",
            "replicate group"
        ]
        rows = []
        for layout in project.layouts:
            for content in layout.layout_contents:
                row = []
                row.append(layout.geid)
                row.append(content.position)
                if content.guide:
                    row.append(content.guide.name)
                else:
                    row.append('')
                row.append(content.sequencing_barcode)
                row.append(content.sequencing_dna_source)
                row.append(content.sequencing_project_id)
                row.append(content.sequencing_sample_name)
                row.append(content.sequencing_library_type)
                if content.clone:
                    row.append(content.clone.cell_line_name)
                    row.append(content.clone.name)
                    row.append(content.clone.cell_pool)
                else:
                    row.append('')
                    row.append('')
                    row.append('')
                row.append(content.content_type)
                row.append(content.is_control)
                row.append(content.replicate_group)
                rows.append(row)
        return headers, rows

    @view_config(route_name="project", renderer="../templates/project.pt")
    def project(self):
        gepid = self.request.matchdict['gepid']
        project = self.dbsession.query(Project).filter(Project.geid == gepid).one()
        plotter = Plotter(self.dbsession, project.geid)
        # Tables
        project_headers, project_rows = self.get_project_table(project)
        target_headers, target_rows = self.get_target_table(project)
        guide_headers, guide_rows = self.get_guide_table(project)
        amplicon_headers, amplicon_rows = self.get_amplicon_table(project)
        layout_headers, layout_rows = self.get_layout_table(project)

        return_map = {'project': project,
                'title': "GenEditID",
                'subtitle': "Project: {}".format(project.geid),
                'info': False,
                'error': False,
                'project_headers': project_headers,
                'project_rows': project_rows,
                'target_headers': target_headers,
                'target_rows': target_rows,
                'guide_headers': guide_headers,
                'guide_rows': guide_rows,
                'amplicon_headers': amplicon_headers,
                'amplicon_rows': amplicon_rows,
                'layout_headers': layout_headers,
                'layout_rows': layout_rows,
                'coverageplot': plotter.coverage_plot(),
                'impactplot': plotter.variant_impact_plot(),
                'heatmapplot': plotter.heatmap_plot(),
            }
        if 'submit_project_data' in self.request.params:
            file_path = None
            try:
                error = False
                file_path = self._upload("layoutfile", ".xlsx")
                self.dbsession.begin_nested()
                loader = ProjectDataLoader(self.dbsession, project.geid, file_path)
                loader.load()
                self.dbsession.commit()
                shutil.copyfile(file_path, os.path.join(project.project_folder, "{}.xlsx".format(project.geid)))
                project_headers, project_rows = self.get_project_table(project)
                target_headers, target_rows = self.get_target_table(project)
                guide_headers, guide_rows = self.get_guide_table(project)
                amplicon_headers, amplicon_rows = self.get_amplicon_table(project)
                layout_headers, layout_rows = self.get_layout_table(project)
                return_map = {'project': project,
                        'title': "GenEditID",
                        'subtitle': "Project: {}".format(project.geid),
                        'info': "Project configuration file has been uploaded successfully.",
                        'error': False,
                        'project_headers': project_headers,
                        'project_rows': project_rows,
                        'target_headers': target_headers,
                        'target_rows': target_rows,
                        'guide_headers': guide_headers,
                        'guide_rows': guide_rows,
                        'amplicon_headers': amplicon_headers,
                        'amplicon_rows': amplicon_rows,
                        'layout_headers': layout_headers,
                        'layout_rows': layout_rows,
                        'coverageplot': plotter.coverage_plot(),
                        'impactplot': plotter.variant_impact_plot(),
                        'heatmapplot': plotter.heatmap_plot(),
                    }
                return HTTPFound(self.request.route_url('project', gepid=project.geid))
            except FinderException as e:
                self.dbsession.rollback()
                error = "Unexpected finder error: {}".format(e)
                self.logger.error(error)
            except LoaderException as e:
                self.dbsession.rollback()
                error = "Unexpected loader error: {}".format(e)
                self.logger.error(error)
            except ValueError as e:
                self.dbsession.rollback()
                error = "Unexpected value error: {}".format(e)
                self.logger.error(error)
            except DBAPIError as e:
                self.dbsession.rollback()
                error = "Unexpected database error: {}".format(e)
                self.logger.error(error)
            except AttributeError as e:
                self.dbsession.rollback()
                error = "Unexpected error: {}".format(e)
                if not file_path:
                    error = "No file selected. Please select a file before uploading."
                self.logger.error(error)
            finally:
                if file_path:
                    try:
                        os.remove(file_path)
                    except OSError:
                        pass
                return_map['error'] = error
                return return_map

        return return_map

    # @view_config(route_name="project_edit", renderer="../templates/project/editproject.pt")
    # def edit_project(self):
    #     id = self.request.matchdict['projectid']
    #     project = self.dbsession.query(Project).filter(Project.id == id).one()
    #     title = "GenEditID"
    #     subtitle = "Project: {}".format(project.geid)
    #     # Project table
    #     project_headers, project_rows = self.get_project_table(project)
    #     comments_error = False
    #     upload_error_project_data = False
    #     upload_error = False
    #     upload_clash = False
    #     plate_headers = [
    #         "layout geid",
    #         "plate geid",
    #         "barcode",
    #         "description",
    #         "data"]
    #     plate_rows = []
    #     plates = self.dbsession.query(Plate).join(Layout).join(Project).filter(Project.geid == project.geid).all()
    #     for plate in plates:
    #         row = []
    #         row.append(plate.layout.geid)
    #         row.append(plate.geid)
    #         row.append(plate.barcode)
    #         row.append(plate.description)
    #         row.append(plate.plate_type)
    #         plate_rows.append(row)
    #     if 'submit_project_data' in self.request.params:
    #         file_path = None
    #         try:
    #             file_path = self._upload("layoutfile", ".xlsx")
    #             self.dbsession.begin_nested()
    #             loader = ProjectDataLoader(self.dbsession, project.geid, file_path)
    #             loader.load()
    #             self.dbsession.commit()
    #             url = self.request.route_url('project_edit', projectid=project.id)
    #             shutil.copyfile(file_path, os.path.join(project.project_folder, "{}.xlsx".format(project.geid)))
    #             return HTTPFound(url)
    #         except LoaderException as e:
    #             self.dbsession.rollback()
    #             self.logger.error("Unexpected loader error: {}".format(e))
    #             upload_error_project_data = str(e)
    #         except ValueError as e:
    #             self.dbsession.rollback()
    #             self.logger.error("Unexpected value error: {}".format(e))
    #             upload_error_project_data = str(e)
    #         except DBAPIError as e:
    #             self.dbsession.rollback()
    #             self.logger.error("Unexpected database error: {}".format(e))
    #             upload_error_project_data = str(e)
    #         except AttributeError as e:
    #             self.dbsession.rollback()
    #             self.logger.error("Unexpected error: {}".format(e))
    #             upload_error_project_data = 'No file selected'
    #         finally:
    #             if file_path:
    #                 try:
    #                     os.remove(file_path)
    #                 except OSError:
    #                     pass
    #     edit_form = self.projects_form("Update")
    #     if 'submit_comments' in self.request.params:
    #         fields = self.request.POST.items()
    #         try:
    #             appstruct = edit_form.validate(fields)
    #         except deform.ValidationFailure as e:
    #             return dict(title=title,
    #                         subtitle=subtitle,
    #                         projectid=project.id,
    #                         project=project,
    #                         project_headers=project_headers,
    #                         project_rows=project_rows,
    #                         comments_error=str(e.error),
    #                         plate_headers=plate_headers,
    #                         plate_rows=plate_rows,
    #                         platelayouts=[p.geid for p in plates],
    #                         platetypes=['abundance (icw)', 'growth (incu)'],
    #                         error_project_data=upload_error_project_data,
    #                         clash=upload_clash,
    #                         error=upload_error)
    #         self.logger.debug("New comments = %s" % appstruct['comments'])
    #         project.comments = appstruct['comments']
    #         url = self.request.route_url('project_edit', projectid=project.id)
    #         return HTTPFound(url)
    #     if 'submit' in self.request.params:
    #         print(self.request.params)
    #         clean_existing = False
    #         try:
    #             clean_existing = self.request.POST['blat']
    #         except KeyError:
    #             pass
    #         file_path = None
    #         try:
    #             file_path = self._upload("datafile")
    #             if file_path and self.request.POST["plate_type"] and self.request.POST["plate_selector"]:
    #                 if self.request.POST["plate_type"].startswith('growth'):
    #                     loader = CellGrowthLoader(self.dbsession, file_path, self.request.POST["plate_selector"])
    #                 if self.request.POST["plate_type"].startswith('abundance'):
    #                     loader = ProteinAbundanceLoader(self.dbsession, file_path, self.request.POST["plate_selector"])
    #                 try:
    #                     loader.load(clean_existing)
    #                     url = self.request.route_url('project_edit', projectid=project.id)
    #                     return HTTPFound(url)
    #                 except ExistingEntityException as e:
    #                     upload_clash = e
    #         except AttributeError as e:
    #             self.dbsession.rollback()
    #             self.logger.error("Unexpected error: {}".format(e))
    #             upload_error = 'No file selected'
    #         except Exception as e:
    #             self.logger.error("Unexpected error: {}".format(e))
    #             upload_error = str(e)
    #         finally:
    #             if file_path:
    #                 try:
    #                     os.remove(file_path)
    #                 except OSError:
    #                     pass
    #     return dict(title=title,
    #                 subtitle=subtitle,
    #                 projectid=project.id,
    #                 project=project,
    #                 project_headers=project_headers,
    #                 project_rows=project_rows,
    #                 comments_error=comments_error,
    #                 plate_headers=plate_headers,
    #                 plate_rows=plate_rows,
    #                 platelayouts=[p.geid for p in plates],
    #                 platetypes=['abundance (icw)', 'growth (incu)'],
    #                 error_project_data=upload_error_project_data,
    #                 clash=upload_clash,
    #                 error=upload_error)

    def _upload(self, property, suffix='.txt'):
        try:
            temp_file_path = ''
            self.logger.debug(self.request.POST[property])
            filename = self.request.POST[property].filename
            filedata = self.request.POST[property].file
            if not filedata:
                return None
            self.logger.debug("Uploaded = %s" % filename)
            file_path = os.path.join('uploads', "{}{}".format(uuid.uuid4(), suffix))
            temp_file_path = file_path + '~'
            filedata.seek(0)
            with open(temp_file_path, 'wb') as output_file:
                shutil.copyfileobj(filedata, output_file)
            os.rename(temp_file_path, file_path)
            statinfo = os.stat(file_path)
            self.logger.info("Uploaded a file of {:d} bytes".format(statinfo.st_size))
        finally:
            try:
                os.remove(temp_file_path)
            except TypeError as e:
                raise(e)
            except OSError:
                pass
        return file_path
