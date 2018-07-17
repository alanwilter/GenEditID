import os
import uuid
import shutil
import logging

import colander
import deform.widget

from collections import OrderedDict

from pyramid.httpexceptions import HTTPFound
from pyramid.view import view_config

from json import JSONEncoder
from sqlalchemy.exc import DBAPIError

from dnascissors.loader import CellGrowthLoader
from dnascissors.loader import ExistingEntityException
from dnascissors.loader import LayoutLoader
from dnascissors.loader import LoaderException
from dnascissors.loader import ProteinAbundanceLoader

from dnascissors.model import ExperimentLayout
from dnascissors.model import Plate
from dnascissors.model import Project
from dnascissors.model import SequencingLibraryContent
from dnascissors.model import Well

from webapp.plots.plotter import Plotter
from webapp.plots.ngsplotter import NGSPlotter
from webapp.clarity import ClaritySubmitter


# See http://docs.pylonsproject.org/projects/pyramid/en/latest/quick_tutorial/forms.html
# File uploads: http://docs.pylonsproject.org/projects/pyramid-cookbook/en/latest/forms/file_uploads.html


class ProjectContent(colander.MappingSchema):
    comments = colander.SchemaNode(colander.String(), title="comments", default='', missing='')


class ProjectViews(object):

    def __init__(self, request):
        self.logger = logging.getLogger(__name__)
        self.request = request
        self.dbsession = request.dbsession
        self.clarity = ClaritySubmitter()

    def projects_form(self, buttonTitle):
        schema = ProjectContent().bind(request=self.request)
        submitButton = deform.form.Button(name='submit_comments', title=buttonTitle)
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
                mismatch_counts = [0] * 4
                mismatch_counts[0] = "{} coding".format(guide.name)
                for mismatch in guide.guide_mismatches:
                    if mismatch.is_off_target_coding_region:
                        if mismatch.number_of_mismatches == 1:
                            mismatch_counts[1] = mismatch.number_of_off_targets
                        elif mismatch.number_of_mismatches == 2:
                            mismatch_counts[2] = mismatch.number_of_off_targets
                        elif mismatch.number_of_mismatches == 3:
                            mismatch_counts[3] = mismatch.number_of_off_targets
                guide_mismatch_rows.append(mismatch_counts)
                mismatch_counts = [0] * 4
                mismatch_counts[0] = "{} non-coding".format(guide.name)
                for mismatch in guide.guide_mismatches:
                    if not mismatch.is_off_target_coding_region:
                        if mismatch.number_of_mismatches == 1:
                            mismatch_counts[1] = mismatch.number_of_off_targets
                        elif mismatch.number_of_mismatches == 2:
                            mismatch_counts[2] = mismatch.number_of_off_targets
                        elif mismatch.number_of_mismatches == 3:
                            mismatch_counts[3] = mismatch.number_of_off_targets
                guide_mismatch_rows.append(mismatch_counts)
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
        sample_data_table_headers = []
        sample_data_table_rows = []
        row_odict = OrderedDict()
        for layout in layouts:
            for well in layout.wells:
                row_odict = OrderedDict()
                row_odict['plate'] = layout.geid
                row_odict['well'] = "{:s}{:02}".format(well.row, well.column)
                row_odict['type'] = 'empty-well'
                row_odict['guides'] = 'no-guide'
                row_odict['zygosity'] = 'empty-well'
                row_odict['protein'] = None
                row_odict['score'] = None
                row_odict['sample'] = None
                row_odict['barcode'] = None
                row_odict['dna source'] = None
                row_odict['consequence'] = None
                row_odict['has off-target'] = None
                row_odict['variant caller'] = None
                row_odict['variant results'] = None
                if well.well_content:
                    row_odict['type'] = well.well_content.content_type
                    row_odict['zygosity'] = 'wt'
                    if well.well_content.guides:
                        row_odict['guides'] = "; ".join(g.name for g in well.well_content.guides)
                for slc in well.sequencing_library_contents:
                    if slc.dna_source == 'fixed cells':
                        row_odict['sample'] = slc.sequencing_sample_name
                        row_odict['barcode'] = slc.sequencing_barcode
                        row_odict['dna source'] = slc.dna_source
                        if slc.mutation_summaries:
                            row_odict['zygosity'] = slc.mutation_summaries[0].zygosity
                            row_odict['consequence'] = slc.mutation_summaries[0].consequence
                            row_odict['has off-target'] = slc.mutation_summaries[0].has_off_target
                            row_odict['variant caller'] = slc.mutation_summaries[0].variant_caller_presence
                            row_odict['score'] = slc.mutation_summaries[0].score
                        if slc.variant_results:
                            row_odict['variant results'] = "; ".join(vr.variant_str_summary for vr in slc.variant_results)
                if well.abundances:
                    if well.abundances[0].ratio_800_700:
                        row_odict['protein'] = "{0:.3f}".format(well.abundances[0].ratio_800_700)
                sample_data_table_rows.append(list(row_odict.values()))
        if row_odict:
            sample_data_table_headers = row_odict.keys()

        return dict(project=project,
                    title="Genome Editing Core",
                    subtitle="Project: {}".format(project.geid),
                    can_sequence=not self.does_pool_exist_in_clarity(project),
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
        comments_error = False
        upload_error_project_data = False
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
        if 'submit_project_data' in self.request.params:
            file_path = None
            try:
                file_path = self._upload("layoutfile", ".xlsx")
                self.dbsession.begin_nested()
                loader = LayoutLoader(self.dbsession, file_path)
                loader.load_project_data(project.id)
                self.dbsession.commit()
                url = self.request.route_url('project_edit', projectid=project.id)
                return HTTPFound(url)
            except LoaderException as e:
                self.dbsession.rollback()
                self.logger.error("Unexpected loader error: {}".format(e))
                upload_error_project_data = str(e)
            except ValueError as e:
                self.dbsession.rollback()
                self.logger.error("Unexpected value error: {}".format(e))
                upload_error_project_data = str(e)
            except DBAPIError as e:
                self.dbsession.rollback()
                self.logger.error("Unexpected database error: {}".format(e))
                upload_error_project_data = str(e)
            except AttributeError as e:
                self.dbsession.rollback()
                self.logger.error("Unexpected error: {}".format(e))
                upload_error_project_data = 'No file selected'
            finally:
                if file_path:
                    try:
                        os.remove(file_path)
                    except OSError:
                        pass
        edit_form = self.projects_form("Update")
        if 'submit_comments' in self.request.params:
            fields = self.request.POST.items()
            try:
                appstruct = edit_form.validate(fields)
            except deform.ValidationFailure as e:
                return dict(title=title,
                            subtitle=subtitle,
                            projectid=project.id,
                            project=project,
                            comments_error=str(e.error),
                            plate_headers=plate_headers,
                            plate_rows=plate_rows,
                            platelayouts=[p.geid for p in plates],
                            platetypes=['abundance (icw)', 'growth (incu)'],
                            error_project_data=upload_error_project_data,
                            clash=upload_clash,
                            error=upload_error)
            self.logger.debug("New comments = %s" % appstruct['comments'])
            project.comments = appstruct['comments']
            url = self.request.route_url('project_edit', projectid=project.id)
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
            except AttributeError as e:
                self.dbsession.rollback()
                self.logger.error("Unexpected error: {}".format(e))
                upload_error = 'No file selected'
            except Exception as e:
                self.logger.error("Unexpected error: {}".format(e))
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
                    comments_error=comments_error,
                    plate_headers=plate_headers,
                    plate_rows=plate_rows,
                    platelayouts=[p.geid for p in plates],
                    platetypes=['abundance (icw)', 'growth (incu)'],
                    error_project_data=upload_error_project_data,
                    clash=upload_clash,
                    error=upload_error)

    @view_config(route_name="project_sequence", renderer="../templates/project/sequenceproject.pt")
    def sequence_project(self):
        id = self.request.matchdict['projectid']
        geproject = self.dbsession.query(Project).filter(Project.id == id).one()

        view_map = dict(title="Genome Editing Core",
                        subtitle="Submit Project {} to Genomics".format(geproject.geid),
                        geproject=geproject,
                        can_sequence=False,
                        error=None,
                        submisson_complete=False)

        slx = self.get_sequencing_id(geproject)
        if not slx:
            view_map['error'] = 'There is no sequencing information for the project.'\
                                'This needs to have been defined in the initial submission.'
            return view_map

        if self.clarity.does_slx_exist(slx):
            view_map['error'] = 'There are already samples in Clarity Genomics LiMS for {}.'.format(slx)
            return view_map

        view_map['can_sequence'] = True

        plate_id = self.request.params.get('plate_id')
        lab_id = self.request.params.get('lab_id')
        researcher_id = self.request.params.get('researcher_id')
        project_source = self.request.params.get('project_source')
        project_id = self.request.params.get('project_id')
        new_project_name = self.request.params.get('new_project_name')

        project_map = self.clarity.get_lab_researcher_project_map()
        json = JSONEncoder(separators=(',', ':')).encode(project_map)

        plates = self.dbsession\
                     .query(Plate)\
                     .join(Plate.experiment_layout)\
                     .join(ExperimentLayout.project)\
                     .filter(Project.id == geproject.id)\
                     .order_by(Plate.geid)\
                     .all()

        sample_count = self.dbsession\
                           .query(SequencingLibraryContent)\
                           .join(SequencingLibraryContent.well)\
                           .join(Well.experiment_layout)\
                           .join(ExperimentLayout.project)\
                           .filter(Project.id == geproject.id)\
                           .count()

        # print("lab_id = {}, researcher_id = {}, project_id = {}".format(lab_id, researcher_id, project_id))

        if 'go_sequence' in self.request.params:
            print("lab_id = {}, researcher_id = {}, submit to = {}, project_id = {}, new project = {}".format(lab_id, researcher_id, project_source, project_id, new_project_name))

            if not plate_id:
                view_map['error'] = "There is plate given to submit from. Go back and select a plate."
                return view_map

            plate = self.dbsession.query(Plate).get(plate_id)

            if not plate:
                view_map['error'] = "There is no such plate in the database. It can only have been removed since the page was loaded."
                return view_map

            if project_source == 'new':
                if not new_project_name:
                    view_map['error'] = "There is no name given for the new project. Go back and supply a project name."
                    return view_map
                if self.clarity.does_project_exist(new_project_name):
                    view_map['error'] = "There is already a project called '{}' in Clarity. Go back and give a different project name."
                    return view_map

                try:
                    clarity_project = self.clarity.create_project(researcher_id, new_project_name)
                    project_id = clarity_project.limsid
                    self.logger.info("Created new project '{}' in Clarity. LIMS id is {}.".format(new_project_name, project_id))
                    project_source = 'existing'
                    new_project_name = None
                except Exception as e:
                    self.logger.error("Error creating new project in Clarity: {}".format(e))
                    view_map['error'] = "There has been a failure creating the new project. {}".format(str(e))
                    return view_map

            self.logger.info("Submitting {} samples into project {} as {} from plate {}.".format(sample_count, project_id, slx, plate.geid))

            try:
                self.clarity.submit_samples(project_id, plate)
                view_map['submisson_complete'] = True
                self.logger.info("Submission complete.")
            except Exception as e:
                self.logger.error("Error submitting samples to Clarity: {}".format(e))
                view_map['error'] = "There has been a failure submitting the samples. {}".format(str(e))
                return view_map

        view_map['jsonmap'] = json
        view_map['plates'] = plates
        view_map['plate_id'] = plate_id
        view_map['lab_id'] = lab_id
        view_map['researcher_id'] = researcher_id
        view_map['project_source'] = project_source
        view_map['project_id'] = project_id
        view_map['new_project_name'] = new_project_name
        view_map['slx'] = slx
        view_map['sample_count'] = sample_count
        view_map['submit_time'] = sample_count / 2
        return view_map

    def get_sequencing_id(self, geproject):
        sequencing_content = self.dbsession\
                                 .query(SequencingLibraryContent)\
                                 .join(SequencingLibraryContent.well)\
                                 .join(Well.experiment_layout)\
                                 .join(ExperimentLayout.project)\
                                 .filter(Project.id == geproject.id)\
                                 .first()
        if sequencing_content:
            return sequencing_content.sequencing_library.slxid
        return None

    def does_pool_exist_in_clarity(self, geproject):
        slx = self.get_sequencing_id(geproject)
        return self.clarity.does_slx_exist(slx)

    def _upload(self, property, suffix='.txt'):
        try:
            temp_file_path = ''
            self.logger.debug(self.request.POST[property])
            filename = self.request.POST[property].filename
            filedata = self.request.POST[property].file
            if not filedata:
                return None
            self.logger.debug("Uploaded = %s" % filename)
            file_path = os.path.join('uploads/', "{}{}".format(uuid.uuid4(), suffix))
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
