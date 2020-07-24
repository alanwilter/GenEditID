import os
import logging
import datetime
import csv
import pandas
import requests
import json

from sqlalchemy.orm.exc import NoResultFound

from geneditid.config import cfg
from geneditid.finder import AmpliconFinder

# reference data
from geneditid.model import Genome
# project data
from geneditid.model import Project
from geneditid.model import Target
from geneditid.model import Guide
from geneditid.model import GuideMismatch
from geneditid.model import Amplicon
from geneditid.model import Primer
from geneditid.model import Clone
from geneditid.model import Layout
from geneditid.model import LayoutContent
from geneditid.model import Plate
# analysis data
from geneditid.model import CellGrowth
from geneditid.model import ProteinAbundance


# --------------------------------------------------------------------------------
# Loader exceptions
# --------------------------------------------------------------------------------
class ExistingEntityException(Exception):

    def __init__(self, objtype, key, msg=None):
        if msg is None:
            msg = "Already have a {} identified by {}.".format(type(objtype).__name__, str(key))
        super(ExistingEntityException, self).__init__(msg)
        self.object_type = type
        self.key = key
        self.message = msg

    def __str__(self):
        return self.message


class LoaderException(Exception):

    def __init__(self, msg=None):
        if msg is None:
            msg = "Loader error."
        super(LoaderException, self).__init__(msg)
        self.message = msg

    def __str__(self):
        return self.message


# --------------------------------------------------------------------------------
# Loader class
# --------------------------------------------------------------------------------
class Loader:

    def get_value(self, text):
        if not text:
            return None
        elif pandas.isnull(text):
            return None
        elif str(text) == 'nan':
            return None
        elif text == '':
            return None
        return text

    def get_string(self, text, max_length=-1):
        if self.get_value(text):
            value = str(text).strip()
            if max_length >= 0 and len(value) > max_length:
                if max_length > 20:
                    return value[:max_length - 3] + "..."
                else:
                    return value[:max_length]
            return value

    def get_int(self, text):
        if self.get_value(text):
            try:
                return int(text)
            except ValueError:
                return None
        return None

    def get_float(self, text):
        if self.get_value(text):
            try:
                return float(text)
            except ValueError:
                return None
        return None

    def get_date(self, text, format='%Y%m%d'):
        if self.get_value(text):
            try:
                return datetime.strptime(str(text), '%Y%m%d')
            except ValueError:
                return None
        return None

    def to_strand(self, value, i):
        if self.get_value(value):
            if value in ['+', 'positive', 'p', 'forward', 'f']:
                return 'forward'
            elif value in ['-', 'negative', 'n', 'reverse', 'r']:
                return 'reverse'
            raise ValueError('Strand must be "forward" or "reverse", not "{}" (row {})'.format(str(value), i))

    def to_dna_feature(self, value, i):
        if self.get_value(value):
            if value in ['gene', 'precursor', 'non-coding']:
                return value
            raise ValueError('DNA feature must be "gene", "precursor", "non-coding", not "{}" (row {})'.format(str(value), i))

    def to_content_type(self, value, i):
        if self.get_value(value):
            value = str(value).lower()
            if value in ['wt', 'wildtype', 'wild-type']:
                return 'wild-type'
            elif value in ['ko', 'knockout', 'knock-out']:
                return 'knock-out'
            elif value in ['bg', 'background']:
                return 'background'
            elif value in ['nm', 'normalisation', 'normaliser']:
                return 'normalisation'
            elif value in ['sm', 'sample']:
                return 'sample'
            elif value == 'empty':
                return 'empty'
            elif value == 'empty-vector':
                return 'empty-vector'
            raise LoaderException("Well content type not recognised: {} (row {})".format(str(value), i))
        else:
            return 'empty'

    def get_gene_id(self, gene_symbol, species='homo_sapiens'):
        return gene_symbol
        # code used to retrive gene_id from ensembl
        # issue with ensembl on 10/07/2020 unable to make it work again
        # species = species.lower().replace(' ', '_')
        # url = ("https://rest.ensembl.org/xrefs/symbol/{}/{}?").format(species, gene_symbol)
        # r = requests.get(url, headers={"Content-Type": "application/json"})
        # gene = json.loads(r.text)
        # if gene:
        #     return gene[0]['id']
        # else:
        #     return gene_symbol


# --------------------------------------------------------------------------------
# RefLoader class
# --------------------------------------------------------------------------------
class RefLoader(Loader):

    GENOMES = cfg['GENOMES']

    def __init__(self, dbsession):
        self.log = logging.getLogger(__name__)
        self.dbsession = dbsession

    def load_genomes(self):
        for genome_name in self.GENOMES:
            species = genome_name.split('[')[0][:-1]
            assembly = genome_name.split('[')[1][:-1]
            # Find or create genome
            genome = self.dbsession.query(Genome).filter(Genome.assembly == assembly).first()
            if not genome:
                genome = Genome(species=species, assembly=assembly)
                self.dbsession.add(genome)
                self.log.info('Created genome {}'.format(genome.assembly))


# --------------------------------------------------------------------------------
# ProjectLoader class
# --------------------------------------------------------------------------------
class ProjectLoader(Loader):

    def __init__(self, dbsession):
        self.log = logging.getLogger(__name__)
        self.dbsession = dbsession
        self.project_geid = None
        self.project = None

    def set_next_project_geid(self):
        last_project = self.dbsession.query(Project).order_by(Project.id.desc()).first()
        if last_project:
            self.project_geid = "GEP{:05d}".format(int(last_project.geid[3:]) + 1)
        else:
            self.project_geid = "GEP00001"

    def create_project(self, name):
        self.set_next_project_geid()
        project = Project()
        project.geid = self.project_geid
        project.name = 'pytest project'
        project.start_date = datetime.date.today()
        self.dbsession.add(project)
        self.project = project

    def delete_project(self, project_geid_todelete):
        self.project_geid = project_geid_todelete
        self.project = self.dbsession.query(Project).filter(Project.geid == project_geid_todelete).first()
        if not self.project:
            raise LoaderException("Project {} not found".format(project_geid_todelete))
        self.log.info("Removing project {} and its associated data.".format(project_geid_todelete))
        self.dbsession.delete(self.project)
        self.log.info('Project {} deleted.'.format(self.project.geid))

    def reset_project(self, project_geid_todelete):
        self.log.info("Reseting project {} and deleting its associated data.".format(project_geid_todelete))
        self.delete_project(project_geid_todelete)
        self.dbsession.add(self.project)
        self.log.info('Project {} named {} has been reset.'.format(self.project.geid, self.project.name))


# --------------------------------------------------------------------------------
# ProjectDataLoader class
# --------------------------------------------------------------------------------
class ProjectDataLoader(Loader):

    def __init__(self, dbsession, project_geid, workbook_file):
        self.log = logging.getLogger(__name__)
        self.dbsession = dbsession
        self.project = self.dbsession.query(Project).filter(Project.geid == project_geid).first()
        if not self.project:
            raise LoaderException("Project {} not found".format(project_geid))
        self.xls = pandas.ExcelFile(workbook_file)
        self.genome = None
        self.finder = AmpliconFinder(self.dbsession, self.project.geid)
        # create project folder
        if not os.path.exists(self.project.project_folder):
            os.makedirs(self.project.project_folder)
        # rename amplicount.csv if exists
        if os.path.exists(os.path.join(self.project.project_folder, 'amplicount.csv')):
            nb_files = len(glob.glob(os.path.join(self.project.project_folder, 'amplicount.csv.*')))
            os.rename(os.path.join(self.project.project_folder, 'amplicount.csv'), os.path.join(self.project.project_folder, 'amplicount.csv.{}'.format(nb_files + 1)))


    def load(self):
        project_loader = ProjectLoader(self.dbsession)
        project_loader.reset_project(self.project.geid)
        self.load_targets()
        self.load_guides()
        self.load_guide_mismatches()
        self.load_amplicon()
        self.load_layout()
        self.load_plates()
        self.finder.write_amplicount_config_file()

    def check_mandatory_fields(self, sheet_name, sheet, mandatory_fields):
        for field in mandatory_fields:
            if field in sheet.columns[sheet.isna().any()].tolist():
                raise LoaderException('Column {} needs a value in {} tab'.format(field, sheet_name))

    def load_targets(self):
        sheet = self.xls.parse('Target')
        if sheet.empty:
            raise LoaderException('Target tab is empty, it must be filled')
        mandatory_fields = ['target_name',
                            'target_genome',
                            'target_gene_id',
                            'target_chrom',
                            'target_start',
                            'target_end',
                            'target_strand']
        self.check_mandatory_fields('Target', sheet, mandatory_fields)
        for i, row in enumerate(sheet.itertuples(), 1):
            self.genome = self.dbsession.query(Genome).filter(Genome.assembly == row.target_genome.split('[')[1][:-1]).first()
            if not self.genome:
                raise LoaderException('Genome {} not found (Target tab, row {})'.format(row.target_genome, i))
            target = Target(project=self.project, genome=self.genome)
            target.name = row.target_name
            target.gene_id = self.get_gene_id(row.target_gene_id, self.genome.species)
            self.log.debug('Target gene_id {}'.format(target.gene_id))
            target.chromosome = str(row.target_chrom)
            target.start = int(row.target_start)
            target.end = int(row.target_end)
            target.strand = self.to_strand(row.target_strand, i)
            target.description = self.get_string(row.target_description, 1024)
            self.dbsession.add(target)
            self.log.info('Target {} created'.format(target.name))

    def load_guides(self):
        sheet = self.xls.parse('Guide')
        if sheet.empty:
            raise LoaderException('Guide tab is empty, it must be filled')
        mandatory_fields = ['target_name',
                            'guide_name',
                            'guide_sequence']
        self.check_mandatory_fields('Guide', sheet, mandatory_fields)
        for i, row in enumerate(sheet.itertuples(), 1):
            target = self.dbsession.query(Target).filter(Target.name == row.target_name).filter(Target.project == self.project).first()
            if not target:
                raise LoaderException('Target {} not found (Guide tab, row {})'.format(row.target_name, i))
            guide = Guide(project=self.project, target=target, genome=self.genome)
            guide.name = row.guide_name
            guide.guide_sequence = row.guide_sequence
            guide.pam_sequence = row.guide_pam_sequence
            if pandas.notnull(row.guide_activity):
                guide.activity = int(row.guide_activity)
            if pandas.notnull(row.guide_exon):
                guide.exon = int(row.guide_exon)
            guide.nuclease = row.guide_nuclease
            self.dbsession.add(guide)
            self.log.info('Guide {} created'.format(guide.name))

    def load_guide_mismatches(self):
        sheet = self.xls.parse('GuideMismatches')
        mandatory_fields = ['guide_name',
                            'is_off_target_coding_region',
                            'number_of_mismatches',
                            'number_of_off_targets']
        if not sheet.empty:
            self.check_mandatory_fields('GuideMismatches', sheet, mandatory_fields)
        for i, row in enumerate(sheet.itertuples(), 1):
            guide = self.dbsession.query(Guide)\
                                  .filter(Guide.project == self.project)\
                                  .filter(Guide.name == row.guide_name)\
                                  .one()
            if not guide:
                raise LoaderException('Guide {} not found (GuideMismatches tab, row {})'.format(row.guide_name, i))
            guide_mismatch = GuideMismatch(guide=guide)
            guide_mismatch.is_off_target_coding_region = bool(row.is_off_target_coding_region)
            guide_mismatch.number_of_mismatches = int(row.number_of_mismatches)
            guide_mismatch.number_of_off_targets = int(row.number_of_off_targets)
            self.dbsession.add(guide_mismatch)
            self.log.info('Guide mismatch ({}, {}, {}) for {} created'.format(guide_mismatch.is_off_target_coding_region, guide_mismatch.number_of_mismatches, guide_mismatch.number_of_off_targets, guide.name))

    def load_amplicon(self):
        sheet = self.xls.parse('Amplicon')
        if sheet.empty:
            raise LoaderException('Amplicon tab is empty, it must be filled')
        mandatory_fields = ['guide_name',
                            'experiment_type',
                            'guide_location',
                            'is_on_target',
                            'dna_feature',
                            'chrom',
                            'forward_primer_sequence',
                            'forward_primer_start',
                            'forward_primer_end',
                            'reverse_primer_sequence',
                            'reverse_primer_start',
                            'reverse_primer_end']
        self.check_mandatory_fields('Amplicon', sheet, mandatory_fields)
        for i, row in enumerate(sheet.itertuples(), 1):
            guide = self.dbsession.query(Guide)\
                                  .filter(Guide.project == self.project)\
                                  .filter(Guide.name == row.guide_name)\
                                  .one()
            if not guide:
                raise LoaderException('Guide {} not found (Amplicon tab, row {})'.format(row.guide_name, i))
            # create the amplicon
            amplicon = Amplicon(project=self.project, genome=self.genome, guide=guide)
            amplicon.dna_feature = self.to_dna_feature(row.dna_feature, i)
            amplicon.chromosome = str(row.chrom)
            amplicon.start = int(row.forward_primer_start)
            amplicon.end = int(row.reverse_primer_end)
            amplicon.experiment_type = row.experiment_type.lower()
            amplicon.guide_location = int(row.guide_location)
            amplicon.is_on_target = bool(row.is_on_target)
            if pandas.notnull(row.score):
                amplicon.score = int(row.score)
            amplicon.description = self.get_string(row.description, 1024)
            self.dbsession.add(amplicon)
            self.log.info('Amplicon {} created'.format(amplicon.name))
            # re-order primer coordinates: reverse primer coordinnates after forward primer
            primer_coords = [int(row.forward_primer_start), int(row.forward_primer_end), int(row.reverse_primer_start), int(row.reverse_primer_end)]
            primer_coords.sort()
            # TODO swap sequences when order different than submitted
            # create forward primer
            forward_primer = Primer(amplicon=amplicon, genome=self.genome)
            forward_primer.sequence = row.forward_primer_sequence.upper()
            forward_primer.strand = 'forward'
            forward_primer.start = primer_coords[0]
            forward_primer.end = primer_coords[1]
            self.dbsession.add(forward_primer)
            self.log.info('Forward primer {} for amplicon {} created'.format(forward_primer.sequence, amplicon.name))
            # create reverse primer
            reverse_primer = Primer(amplicon=amplicon, genome=self.genome)
            reverse_primer.sequence = row.reverse_primer_sequence.upper()
            reverse_primer.strand = 'reverse'
            reverse_primer.start = primer_coords[2]
            reverse_primer.end = primer_coords[3]
            self.dbsession.add(reverse_primer)
            self.log.info('Reverse primer {} for amplicon {} created'.format(reverse_primer.sequence, amplicon.name))

    def load_layout(self):
        sheet = self.xls.parse('Layout')
        if sheet.empty:
            raise LoaderException('Layout tab is empty, it must be filled')
        mandatory_fields = ['layout_id',
                            'well_position']
        self.check_mandatory_fields('Layout', sheet, mandatory_fields)
        for i, row in enumerate(sheet.itertuples(), 1):
            # Find or create layout
            layout = self.dbsession.query(Layout)\
                                   .filter(Layout.geid == row.layout_id)\
                                   .filter(Layout.project == self.project)\
                                   .first()
            if not layout:
                layout = Layout(project=self.project)
                layout.geid = row.layout_id
                self.dbsession.add(layout)
                self.log.info('Layout {} in project {} created.'.format(layout.geid, self.project.geid))
            # Empty layout content when no sequencing_barcode
            if not self.get_value(row.sequencing_barcode):
                layout_content = LayoutContent(layout=layout)
                layout_content.row = row.well_position[0]
                layout_content.column = int(row.well_position[1:])
                layout_content.content_type = 'empty'
                self.dbsession.add(layout_content)
                self.log.info('Empty LayoutContent in position {}{} in layout {} created'.format(layout_content.row, layout_content.column, layout_content.layout.geid))
            # Layout content
            else:
                if not self.get_value(row.sequencing_barcode):
                    raise LoaderException('sequencing_barcode need a value in Layout tab, row {} when not empty'.format(i))
                # Clone
                clone = None
                if self.get_value(row.clone_name):
                    clone = self.dbsession.query(Clone).filter(Clone.name == str(row.clone_name))\
                                                       .filter(Clone.cell_pool == self.get_string(row.cell_pool))\
                                                       .filter(Clone.cell_line_name == str(row.cell_line_name))\
                                                       .filter(Clone.project == self.project)\
                                                       .first()
                    if not clone:
                        clone = Clone(project=self.project)
                        clone.cell_line_name = str(row.cell_line_name)
                        clone.name = str(row.clone_name)
                        clone.cell_pool=self.get_string(row.cell_pool)
                        self.dbsession.add(clone)
                        self.log.info('Clone {} in cell line {} pool {}'.format(clone.name, clone.cell_line_name, clone.cell_pool))
                # Layout content
                layout_content = LayoutContent(layout=layout)
                layout_content.clone = clone
                layout_content.row = row.well_position[0]
                layout_content.column = int(row.well_position[1:])
                layout_content.sequencing_project_id = self.get_value(row.sequencing_project_id)
                layout_content.sequencing_library_type = self.get_value(row.sequencing_library_type)
                layout_content.sequencing_barcode = self.get_value(row.sequencing_barcode)
                layout_content.sequencing_dna_source = self.get_value(row.sequencing_dna_source)
                layout_content.sequencing_sample_name = self.get_value(row.sequencing_sample_name)
                layout_content.content_type = self.get_value(row.content_type)
                layout_content.is_control = self.get_value(row.is_control)
                layout_content.replicate_group = self.get_value(row.replicate_group)
                self.dbsession.add(layout_content)
                self.log.info('LayoutContent in position {}{} in layout {} created'.format(layout_content.row, layout_content.column, layout_content.layout.geid))

    def load_plates(self):
        sheet = self.xls.parse('Plate')
        mandatory_fields = ['layout_id',
                            'plate_name']
        if not sheet.empty:
            self.check_mandatory_fields('Plate', sheet, mandatory_fields)
        for i, row in enumerate(sheet.itertuples(), 1):
            layout = self.dbsession.query(Layout)\
                                   .filter(Layout.geid == row.layout_id)\
                                   .filter(Layout.project == self.project)\
                                   .first()
            if not layout:
                raise LoaderException('Layout {} not found (Plate tab, row {})'.format(row.layout_id, i))
            check_plate = self.dbsession.query(Plate)\
                                        .filter(Plate.layout == layout)\
                                        .filter(Plate.name == row.plate_name)\
                                        .first()
            if check_plate:
                raise LoaderException('Plate {} already exists for this layout {}, please assign a unique value (Plate tab, row {})'.format(row.plate_name, row.layout_id, i))
            plate = Plate(layout=layout)
            plate.name = row.plate_name
            plate_barcode = self.dbsession.query(Plate).filter(Plate.barcode == row.plate_barcode).first()
            if plate_barcode:
                raise LoaderException('Plate barcode {} already exists, please assign a unique value (Plate tab, row {})'.format(row.plate_barcode, i))
            plate.barcode = row.plate_barcode
            plate.description = row.plate_description
            self.dbsession.add(plate)
            self.log.info('Plate {} in layout {} created'.format(plate.name, plate.layout.geid))


# --------------------------------------------------------------------------------
# ProteinAbundanceLoader class
# --------------------------------------------------------------------------------
class ProteinAbundanceLoader(Loader):

    def __init__(self, dbsession, csv_file, plate_id):
        self.log = logging.getLogger(__name__)
        self.dbsession = dbsession
        self.csv_file = csv_file
        self.plate_id = plate_id
        self.plate = self.dbsession.query(Plate).filter(Plate.geid == plate_id).one()
        if not self.plate:
            raise LoaderException('No plate {}'.format(plate_id))

    def clean(self):
        self.dbsession.query(ProteinAbundance).filter(ProteinAbundance.plate_id == self.plate.id).delete()
        self.dbsession.flush()

    def load(self, clean_if_exists=False):
        if self.plate.is_abundance_plate:
            if clean_if_exists:
                self.clean()
            else:
                raise ExistingEntityException(Plate, self.plate.geid, "Already have protein abundance data for this plate {}. Will not overwrite it.".format(self.plate.geid))
        self.log.info("Loading InCell Western signal for plate {} from {}".format(self.plate_id, self.csv_file))
        with open(self.csv_file, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            channel = None
            row = None
            a = ord('A')
            for line in reader:
                if line[0] == '':
                    channel = None
                    row = None
                elif line[0] == 'Total Channel 700':
                    channel = 700
                    row = 0
                elif line[0] == 'Total Channel 800':
                    channel = 800
                    row = 0
                elif channel:
                    rowchar = chr(a + row)
                    for column in range(1, 13):
                        well = self.dbsession.query(Well)\
                                           .filter(Well.layout == self.plate.layout)\
                                           .filter(Well.row == rowchar)\
                                           .filter(Well.column == column)\
                                           .first()
                        if not well:
                            raise LoaderException("There is no well at {}{} on plate {}".format(rowchar, column, self.plate.geid))
                        self.log.info("Setting well {}{} on {} for ICW {}".format(well.row, well.column, self.plate.geid, channel))

                        abundance = self.dbsession.query(ProteinAbundance)\
                                                .filter(ProteinAbundance.well == well)\
                                                .filter(ProteinAbundance.plate == self.plate)\
                                                .first()
                        if not abundance:
                            abundance = ProteinAbundance(well=well, plate=self.plate)
                            self.dbsession.add(abundance)
                        if channel == 700:
                            abundance.intensity_channel_700 = self.get_float(line[column - 1])
                        elif channel == 800:
                            abundance.intensity_channel_800 = self.get_float(line[column - 1])
                    row = row + 1


# --------------------------------------------------------------------------------
# CellGrowthLoader class
# --------------------------------------------------------------------------------
class CellGrowthLoader(Loader):

    def __init__(self, dbsession, csv_file, plate_id):
        self.log = logging.getLogger(__name__)
        self.dbsession = dbsession
        self.csv_file = csv_file
        self.plate_id = plate_id
        self.plate = self.dbsession.query(Plate).filter(Plate.geid == plate_id).one()
        if not self.plate:
            raise LoaderException('No plate {}'.format(plate_id))

    def clean(self):
        self.dbsession.query(CellGrowth).filter(CellGrowth.plate_id == self.plate.id).delete()
        self.dbsession.flush()

    def load(self, clean_if_exists=False):
        if self.plate.is_growth_plate:
            if clean_if_exists:
                self.clean()
            else:
                raise ExistingEntityException(Plate, self.plate.geid, "Already have cell growth data for this plate {}. Will not overwrite it.".format(self.plate.geid))

        self.log.info("Loading Incucyte growth information for plate {} from {}".format(self.plate_id, self.csv_file))

        with open(self.csv_file, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            start = False
            header = None
            for line in reader:
                if len(line) == 0:
                    continue
                if line[0] == 'Date Time':
                    print("START")
                    start = True
                    header = line
                elif start:
                    time = datetime.strptime(line[0], '%d/%m/%Y %H:%M:%S')
                    hour = float(line[1])
                    for column in range(2, len(line)):
                        location = header[column]
                        wellrow = location[0]
                        wellcolumn = int(location[1:])
                        well = self.dbsession.query(Well)\
                                           .filter(Well.layout == self.plate.layout)\
                                           .filter(Well.row == wellrow)\
                                           .filter(Well.column == wellcolumn)\
                                           .first()
                        if not well:
                            raise LoaderException("There is no well at {}{} on plate {}".format(wellrow, wellcolumn, self.plate.geid))
                        value = float(line[column])
                        self.log.info("Setting well {}{} on {} for Incucyte {:f} at {}".format(well.row, well.column, self.plate.geid, value, line[0]))
                        growth = self.dbsession.query(CellGrowth)\
                                             .filter(CellGrowth.plate == self.plate)\
                                             .filter(CellGrowth.well == well)\
                                             .filter(CellGrowth.hours == hour)\
                                             .first()
                        if not growth:
                            growth = CellGrowth(well=well, plate=self.plate)
                            growth.hours = hour
                            growth.timestamp = time
                            self.dbsession.add(growth)
                        growth.confluence_percentage = value
