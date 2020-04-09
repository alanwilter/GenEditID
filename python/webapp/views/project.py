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
        project_headers = [
            "geid",
            "name",
            "type",
            "scientist",
            "group",
            "date",
            "description",
            "comments",
            "abundance data",
            "growth data",
            "ngs data"]
        project_rows = [[project.geid,
                        project.name,
                        project.project_type,
                        project.scientist,
                        project.group,
                        project.start_date,
                        project.description,
                        project.comments,
                        project.is_abundance_data_available,
                        project.is_growth_data_available,
                        project.is_variant_data_available]]
        return project_headers, project_rows

    def get_variant_data_table(self, layouts, variant_caller='VarDict'):
        variant_data_table_headers = []
        variant_data_table_rows = []
        vrow_odict = OrderedDict()
        for layout in layouts:
            for well in layout.wells:
                genome = None
                if well.well_content and len(well.well_content.guides) > 0:
                    genome = well.well_content.guides[0].genome.assembly
                for slc in well.sequencing_library_contents:
                    if slc.variant_results:
                        for v in slc.variant_results:
                            if v.variant_caller == variant_caller:
                                # https://software.broadinstitute.org/software/igv/ControlIGV
                                igv_url = None
                                bam = "http://bioinf-ge001.cri.camres.org/data/{}/idbam/{}.bam".format(layout.project.geid, slc.sequencing_barcode)
                                if genome:
                                    locus = "{}:{}".format(v.chromosome, v.position)
                                    igv_url = "http://localhost:60151/load?file={}&locus={}&genome={}&name={}".format(bam, locus, genome, quote(slc.sequencing_sample_name))

                                vrow_odict['plate'] = layout.geid
                                vrow_odict['well'] = "{:s}{:02}".format(well.row, well.column)
                                vrow_odict['clone'] = well.well_content.clone.name if well.well_content else None
                                vrow_odict['sample'] = slc.sequencing_sample_name
                                vrow_odict['barcode'] = slc.sequencing_barcode
                                vrow_odict['variant_caller'] = v.variant_caller
                                vrow_odict['variant_type'] = v.variant_type
                                vrow_odict['IGV link'] = igv_url
                                vrow_odict['BAM file'] = bam
                                vrow_odict['consequence'] = v.consequence
                                vrow_odict['gene_id'] = v.gene_id
                                vrow_odict['gene'] = v.gene
                                vrow_odict['cdna_effect'] = v.cdna_effect
                                vrow_odict['protein_effect'] = v.protein_effect
                                vrow_odict['codons'] = v.codons
                                vrow_odict['chromosome'] = v.chromosome
                                vrow_odict['position'] = v.position
                                vrow_odict['ref'] = v.ref
                                vrow_odict['alt'] = v.alt
                                vrow_odict['allele_fraction'] = v.allele_fraction
                                vrow_odict['depth'] = v.depth
                                vrow_odict['quality'] = v.quality
                                vrow_odict['amplicon'] = v.amplicon
                                vrow_odict['exon'] = v.exon
                                vrow_odict['intron'] = v.intron
                                vrow_odict['existing_variation'] = v.existing_variation
                                vrow_odict['sift'] = v.sift
                                vrow_odict['polyphen'] = v.polyphen
                                vrow_odict['clinical_significance'] = v.clinical_significance
                                vrow_odict['gmaf'] = v.gmaf
                                vrow_odict['offset_from_primer_end'] = v.offset_from_primer_end
                                vrow_odict['indel_length'] = v.indel_length
                                vrow_odict['forward_context'] = v.forward_context
                                vrow_odict['alleles'] = v.alleles
                                vrow_odict['reverse_context'] = v.reverse_context
                                variant_data_table_rows.append(list(vrow_odict.values()))
        if vrow_odict:
            variant_data_table_headers = vrow_odict.keys()
        return variant_data_table_headers, variant_data_table_rows

    @view_config(route_name="project_view", renderer="../templates/project/viewproject.pt")
    def view_project(self):
        id = self.request.matchdict['projectid']
        project = self.dbsession.query(Project).filter(Project.id == id).one()
        plotter = Plotter(self.dbsession, project.geid)
        ngsplotter = NGSPlotter(self.dbsession, project.geid)
        # Project table
        project_headers, project_rows = self.get_project_table(project)
        # Target table
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
        guide_rows = []
        # Guide mismatch table
        guide_mismatch_headers = [
            "genome region",
            "1",
            "2",
            "3"
        ]
        guide_mismatch_rows = []
        targets = project.targets
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
        # Scores: sample data table
        layouts = project.experiment_layouts
        sample_data_table_headers = []
        sample_data_table_rows = []
        row_odict = OrderedDict()
        for layout in layouts:
            for well in layout.wells:
                row_odict = OrderedDict()
                row_odict['plate'] = layout.geid
                row_odict['well'] = "{:s}{:02}".format(well.row, well.column)
                row_odict['type'] = None
                row_odict['clone'] = None
                row_odict['guides'] = None
                row_odict['slxid'] = None
                row_odict['sample'] = None
                row_odict['barcode'] = None
                row_odict['dna source'] = None
                row_odict['protein'] = None
                row_odict['zygosity'] = None
                row_odict['consequence'] = None
                row_odict['has off-target'] = None
                row_odict['variant caller presence'] = None
                row_odict['score'] = None
                if well.well_content:
                    row_odict['type'] = well.well_content.content_type
                    row_odict['clone'] = well.well_content.clone.name
                    if well.well_content.guides:
                        row_odict['guides'] = "; ".join(g.name for g in well.well_content.guides)
                for slc in well.sequencing_library_contents:
                    row_odict['slxid'] = slc.sequencing_library.slxid
                    row_odict['sample'] = slc.sequencing_sample_name
                    row_odict['barcode'] = slc.sequencing_barcode
                    row_odict['dna source'] = slc.dna_source
                    if slc.mutation_summaries:
                        row_odict['zygosity'] = slc.mutation_summaries[0].zygosity
                        row_odict['consequence'] = slc.mutation_summaries[0].consequence
                        row_odict['has off-target'] = slc.mutation_summaries[0].has_off_target
                        row_odict['variant caller presence'] = slc.mutation_summaries[0].variant_caller_presence
                        row_odict['score'] = slc.mutation_summaries[0].score
                if well.abundances:
                    if well.abundances[0].ratio_800_700:
                        row_odict['protein'] = "{0:.3f}".format(well.abundances[0].ratio_800_700)
                sample_data_table_rows.append(list(row_odict.values()))
        if row_odict:
            sample_data_table_headers = row_odict.keys()
        vvariant_data_table_headers, vvariant_data_table_rows = self.get_variant_data_table(layouts, 'VarDict')
        hvariant_data_table_headers, hvariant_data_table_rows = self.get_variant_data_table(layouts, 'HaplotypeCaller')

        return dict(project=project,
                    title="GenEditID",
                    subtitle="Project: {}".format(project.geid),
                    #cellgrowthplot=plotter.growth_plot(),
                    coverageplot=plotter.coverage_plot(),
                    impactplot=plotter.impact_plot(),
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
                    vvariant_data_table_headers=vvariant_data_table_headers,
                    vvariant_data_table_rows=vvariant_data_table_rows,
                    hvariant_data_table_headers=hvariant_data_table_headers,
                    hvariant_data_table_rows=hvariant_data_table_rows,
                    sample_data_table_headers=sample_data_table_headers,
                    sample_data_table_rows=sample_data_table_rows)

    @view_config(route_name="project_edit", renderer="../templates/project/editproject.pt")
    def edit_project(self):
        id = self.request.matchdict['projectid']
        project = self.dbsession.query(Project).filter(Project.id == id).one()
        title = "GenEditID"
        subtitle = "Project: {}".format(project.geid)
        # Project table
        project_headers, project_rows = self.get_project_table(project)
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
                            project_headers=project_headers,
                            project_rows=project_rows,
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
                    project_headers=project_headers,
                    project_rows=project_rows,
                    comments_error=comments_error,
                    plate_headers=plate_headers,
                    plate_rows=plate_rows,
                    platelayouts=[p.geid for p in plates],
                    platetypes=['abundance (icw)', 'growth (incu)'],
                    error_project_data=upload_error_project_data,
                    clash=upload_clash,
                    error=upload_error)

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
