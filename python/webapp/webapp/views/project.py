import os
import uuid
import shutil
import logging

import colander
import deform.widget

from pyramid.httpexceptions import HTTPFound
from pyramid.view import view_config

from dnascissors.model import Project
from dnascissors.model import ExperimentLayout
from dnascissors.model import Plate

from dnascissors.loader import ExistingEntityException

from dnascissors.loader import ProteinAbundanceLoader
from dnascissors.loader import CellGrowthLoader

from webapp.plots.plotter import Plotter
from webapp.plots.ngsplotter import NGSPlotter


# See http://docs.pylonsproject.org/projects/pyramid/en/latest/quick_tutorial/forms.html
# File uploads: http://docs.pylonsproject.org/projects/pyramid-cookbook/en/latest/forms/file_uploads.html

class ProjectContent(colander.MappingSchema):
    comments = colander.SchemaNode(colander.String(), title="comments")


class ProjectViews(object):

    def __init__(self, request):
        self.logger = logging.getLogger('webapp')
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
        plotter = Plotter(self.dbsession, project.geid)
        ngsplotter = NGSPlotter(self.dbsession, project.geid)
        # Project table
        project_headers = [
            "geid",
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
            "ngs data"]
        project_rows = [[project.geid,
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
                        project.is_variant_data_available]]
        # Target table
        targets = project.targets
        target_headers = [
            "name",
            "species",
            "assembly",
            "gene",
            "chromosome",
            "start",
            "end",
            "strand",
            "description"]
        target_rows = []
        guide_rows = []
        guide_mismatch_rows = []
        for target in targets:
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
            target_rows.append(row)
            for guide in target.guides:
                guide_row = []
                guide_row.append(target.name)
                guide_row.append(guide.genome.species)
                guide_row.append(guide.genome.assembly)
                guide_row.append(guide.name)
                guide_row.append(guide.guide_sequence)
                guide_row.append(guide.pam_sequence)
                guide_row.append(guide.activity)
                guide_row.append(guide.exon)
                guide_row.append(guide.nuclease)
                guide_rows.append(guide_row)
                mismatch_dict = {}
                mismatch_dict['genome region'] = "{} coding".format(guide.name)
                for mismatch in guide.guide_mismatches:
                    if mismatch.is_off_target_coding_region:
                        if mismatch.number_of_mismatches == 1:
                            mismatch_dict['1'] = mismatch.number_of_off_targets
                        elif mismatch.number_of_mismatches == 2:
                            mismatch_dict['2'] = mismatch.number_of_off_targets
                        elif mismatch.number_of_mismatches == 3:
                            mismatch_dict['3'] = mismatch.number_of_off_targets
                guide_mismatch_rows.append(mismatch_dict.values())
                mismatch_dict = {}
                mismatch_dict['genome region'] = "{} non-coding".format(guide.name)
                for mismatch in guide.guide_mismatches:
                    if not mismatch.is_off_target_coding_region:
                        if mismatch.number_of_mismatches == 1:
                            mismatch_dict['1'] = mismatch.number_of_off_targets
                        elif mismatch.number_of_mismatches == 2:
                            mismatch_dict['2'] = mismatch.number_of_off_targets
                        elif mismatch.number_of_mismatches == 3:
                            mismatch_dict['3'] = mismatch.number_of_off_targets
                guide_mismatch_rows.append(mismatch_dict.values())
        # Guide table
        guide_headers = [
            "target name",
            "species",
            "assembly",
            "guide name",
            "guide sequence",
            "pam sequence",
            "activity",
            "exon",
            "nuclease"
        ]
        # Guide mismatch table
        guide_mismatch_headers = [
            "genome region",
            "1",
            "2",
            "3"
        ]
        # Sample analysis: data table
        layouts = project.experiment_layouts
        sample_data_table_headers = [
            "plate",
            "well",
            "guides",
            "type",
            "sample",
            "barcode",
            "zygosity",
            "consequence",
            "has off-target",
            "variant caller",
            "score",
            #"variant types",
            #"allele fractions",
            "protein",
            ]
        sample_data_table_rows = []
        for layout in layouts:
            for well in layout.wells:
                row = []
                row.append(layout.geid)
                row.append("{:s}{:02}".format(well.row, well.column))
                guide_names = 'no-guide'
                content_type = 'empty-well'
                zygosity = 'empty-well'
                if well.well_content:
                    content_type = well.well_content.content_type
                    zygosity = 'wt'
                    if well.well_content.guides:
                        guide_names = "; ".join(g.name for g in well.well_content.guides)
                row.append(guide_names)
                row.append(content_type)
                slc_fixed_cells_found = False
                for slc in well.sequencing_library_contents:
                    if slc.dna_source == 'fixed cells':
                        slc_fixed_cells_found = True
                        row.append(slc.sequencing_sample_name)
                        row.append(slc.sequencing_barcode)
                        if slc.mutation_summaries:
                            row.append(slc.mutation_summaries[0].zygosity)
                            row.append(slc.mutation_summaries[0].consequence)
                            row.append(slc.mutation_summaries[0].has_off_target)
                            row.append(slc.mutation_summaries[0].variant_caller_presence)
                            row.append(slc.mutation_summaries[0].score)
                        else:
                            row.append(zygosity)
                            row.append(None)
                            row.append(None)
                            row.append(None)
                            row.append(None)
                if not slc_fixed_cells_found:
                    row.append(zygosity)
                    row.append(None)
                    row.append(None)
                    row.append(None)
                    row.append(None)
                    row.append(None)
                    row.append(None)
                    # if slc.variant_results:
                    #     row.append("; ".join(vr.variant_type for vr in slc.variant_results))
                    #     row.append("; ".join("{0:.3f}".format(vr.allele_fraction) for vr in slc.variant_results))
                    # else:
                    #     row.append(None)
                    #     row.append(None)
                well_abundance_ratio = None
                if well.abundances:
                    well_abundance_ratio = "{0:.3f}".format(well.abundances[0].ratio_800_700)
                row.append(well_abundance_ratio)
                sample_data_table_rows.append(row)

        return dict(project=project,
                    title="Genome Editing Core",
                    subtitle="Project: {}".format(project.geid),
                    cellgrowthplot=plotter.growth_plot(),
                    proteinabundanceplot=plotter.abundance_plot(),
                    ngsplot=ngsplotter.combined_ngs_plot(),
                    platescoringplot=plotter.plate_scoring_plots(),
                    project_headers=project_headers,
                    project_rows=project_rows,
                    target_headers=target_headers,
                    target_rows=target_rows,
                    guide_headers=guide_headers,
                    guide_rows=guide_rows,
                    guide_mismatch_headers=guide_mismatch_headers,
                    guide_mismatch_rows=guide_mismatch_rows,
                    sample_data_table_headers=sample_data_table_headers,
                    sample_data_table_rows=sample_data_table_rows)

    @view_config(route_name="project_edit", renderer="../templates/project/editproject.pt")
    def edit_project(self):
        upload_error = False
        upload_clash = False
        plate_headers = [
            "layout geid",
            "plate geid",
            "barcode",
            "description",
            "data"]
        plate_rows = []
        id = self.request.matchdict['projectid']
        project = self.dbsession.query(Project).filter(Project.id == id).one()
        title = "Genome Editing Core"
        subtitle = "Project: {}".format(project.geid)
        plates = self.dbsession.query(Plate).join(ExperimentLayout).join(Project).filter(Project.geid == project.geid).all()
        for plate in plates:
            row = []
            row.append(plate.experiment_layout.geid)
            row.append(plate.geid)
            row.append(plate.barcode)
            row.append(plate.description)
            row.append(plate.plate_type)
            plate_rows.append(row)
        edit_form = self.projects_form("Update")
        if 'submit_comments' in self.request.params:
            fields = self.request.POST.items()
            try:
                appstruct = edit_form.validate(fields)
            except deform.ValidationFailure as e:
                return dict(project=project, form=e.render(), title=title, subtitle=subtitle)
            print("New comments = %s" % appstruct['comments'])
            project.comments = appstruct['comments']
            url = self.request.route_url('project_view', projectid=project.id)
            return HTTPFound(url)
        if 'submit' in self.request.params:
            print(self.request.params)
            clean_existing = False
            try:
                clean_existing = self.request.POST['blat']
            except KeyError:
                pass
            file_path = None
            try:
                file_path = self._upload("datafile")
                if file_path and self.request.POST["plate_type"] and self.request.POST["plate_selector"]:
                    if self.request.POST["plate_type"].startswith('growth'):
                        loader = CellGrowthLoader(self.dbsession, file_path, self.request.POST["plate_selector"])
                    if self.request.POST["plate_type"].startswith('abundance'):
                        loader = ProteinAbundanceLoader(self.dbsession, file_path, self.request.POST["plate_selector"])
                    try:
                        loader.load(clean_existing)
                        url = self.request.route_url('project_edit', projectid=project.id)
                        return HTTPFound(url)
                    except ExistingEntityException as e:
                        upload_clash = e
            except Exception as e:
                self.logger.error("Have an unexpected error while uploading data: {}".format(e))
                upload_error = str(e)
            finally:
                if file_path:
                    try:
                        os.remove(file_path)
                    except OSError:
                        pass
        return dict(title=title,
                    subtitle=subtitle,
                    projectid=project.id,
                    project=project,
                    plate_headers=plate_headers,
                    plate_rows=plate_rows,
                    platelayouts=[p.geid for p in plates],
                    platetypes=['abundance (icw)', 'growth (incu)'],
                    clash=upload_clash,
                    error=upload_error)

    def _upload(self, property):
        filename = self.request.POST[property].filename
        filedata = self.request.POST[property].file
        if not filedata:
            return None
        self.logger.debug("Uploaded = %s" % filename)
        file_path = os.path.join('uploads/', "{}{}".format(uuid.uuid4(), '.txt'))
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
