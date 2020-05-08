import logging
import datetime
import csv
import pandas
import requests
import json

from sqlalchemy.orm.exc import NoResultFound

from geneditid.config import cfg

# reference data
from geneditid.model import Genome
from geneditid.model import CellLine
# project data
from geneditid.model import Project
from geneditid.model import Target
from geneditid.model import Guide
from geneditid.model import GuideMismatch
from geneditid.model import Amplicon
from geneditid.model import AmpliconSelection
from geneditid.model import Primer
from geneditid.model import Clone
from geneditid.model import ExperimentLayout
from geneditid.model import Plate
from geneditid.model import WellContent
from geneditid.model import Well
from geneditid.model import SequencingLibrary
from geneditid.model import SequencingLibraryContent
# analysis data
from geneditid.model import CellGrowth
from geneditid.model import ProteinAbundance
from geneditid.model import VariantResult


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
        species = species.lower().replace(' ', '_')
        url = ("https://rest.ensembl.org/xrefs/symbol/{}/{}?").format(species, gene_symbol)
        r = requests.get(url, headers={"Content-Type": "application/json"})
        gene = json.loads(r.text)
        if gene:
            return gene[0]['id']
        else:
            return gene_symbol


# --------------------------------------------------------------------------------
# RefLoader class
# --------------------------------------------------------------------------------
class RefLoader(Loader):

    GENOMES = cfg['GENOMES']

    CELL_LINES = cfg['CELL_LINES']

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

    def load_celllines(self):
        for cell_line_name in self.CELL_LINES:
            # Find or create cell line
            cell_line = self.dbsession.query(CellLine).filter(CellLine.name == cell_line_name).first()
            if not cell_line:
                cell_line = CellLine(name=cell_line_name)
                self.dbsession.add(cell_line)
                self.log.info('Created cell line {}'.format(cell_line.name))


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

    def create_project(self, name, project_type):
        self.set_next_project_geid()
        project = Project()
        project.geid = self.project_geid
        project.name = 'pytest project'
        project.project_type = 'knock-out'
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

    def load(self):
        project_loader = ProjectLoader(self.dbsession)
        project_loader.reset_project(self.project.geid)
        self.load_targets()
        self.load_guides()
        self.load_guide_mismatches()
        self.load_amplicon()
        #self.load_layout()
        self.load_plates()

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
            guide = Guide(target=target, genome=self.genome)
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
                                  .join(Target)\
                                  .join(Project)\
                                  .filter(Project.geid == self.project.geid)\
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
                            'guide_strand',
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
                                  .join(Target)\
                                  .join(Project)\
                                  .filter(Project.geid == self.project.geid)\
                                  .filter(Guide.name == row.guide_name)\
                                  .one()
            if not guide:
                raise LoaderException('Guide {} not found (Amplicon tab, row {})'.format(row.guide_name, i))
            # create the amplicon
            amplicon = Amplicon(project=self.project, genome=self.genome)
            amplicon.dna_feature = self.to_dna_feature(row.dna_feature, i)
            amplicon.chromosome = str(row.chrom)
            amplicon.start = int(row.forward_primer_start)
            amplicon.end = int(row.reverse_primer_end)
            self.dbsession.add(amplicon)
            self.log.info('Amplicon {} created'.format(amplicon.get_name))
            # create link between amplicon and guide
            selection = AmpliconSelection(guide=guide, amplicon=amplicon)
            selection.experiment_type = row.experiment_type
            selection.guide_location = int(row.guide_location)
            selection.guide_strand = self.to_strand(row.guide_strand, i)
            selection.is_on_target = bool(row.is_on_target)
            if pandas.notnull(row.score):
                selection.score = int(row.score)
            selection.description = self.get_string(row.description, 1024)
            self.dbsession.add(selection)
            self.log.info('Amplicon link between amplicon {} and guide {} at {} created'.format(amplicon.get_name, guide.name, selection.guide_location))
            # check reverse primer coordinnates after forward primer
            if not (int(row.forward_primer_start) < int(row.forward_primer_end) < int(row.reverse_primer_start) < int(row.reverse_primer_end)):
                raise LoaderException('Forward primer coordinates must be before reverse ones (Amplicon tab, row {})'.format(i))
            # create forward primer
            forward_primer = Primer(amplicon=amplicon, genome=self.genome)
            forward_primer.sequence = row.forward_primer_sequence
            forward_primer.strand = 'forward'
            forward_primer.start = int(row.forward_primer_start)
            forward_primer.end = int(row.forward_primer_end)
            self.dbsession.add(forward_primer)
            self.log.info('Forward primer {} for amplicon {} created'.format(forward_primer.sequence, amplicon.get_name))
            # create reverse primer
            reverse_primer = Primer(amplicon=amplicon, genome=self.genome)
            reverse_primer.sequence = row.reverse_primer_sequence
            reverse_primer.strand = 'reverse'
            reverse_primer.start = int(row.reverse_primer_start)
            reverse_primer.end = int(row.reverse_primer_end)
            self.dbsession.add(reverse_primer)
            self.log.info('Reverse primer {} for amplicon {} created'.format(reverse_primer.sequence, amplicon.get_name))

    def load_layout(self):
        sheet = self.xls.parse('Layout')
        for i, row in enumerate(sheet.itertuples(), 1):
            guide = None
            clone = None
            content = None
            # May also need to find a guide.
            if self.get_value(row.guide_name):
                try:
                    guide = self.dbsession.query(Guide).\
                                         join(Target).\
                                         join(Project).\
                                         filter(Project.geid == self.project.geid).\
                                         filter(Guide.name == row.guide_name).\
                                         one()
                except NoResultFound:
                    raise LoaderException('Guide "{}" does not exist (row {})'.format(row.guide_name, i))
            # Find or create experiment layout
            if not row.geid:
                raise LoaderException('Experiment layout GEID is required on row {}'.format(i))
            #if not row.geid.split('_')[0] == self.project.geid:
            #    raise LoaderException('Experiment layout GEID {} does not start with project GEID {}'.format(row.geid, self.project.geid))
            layout = self.dbsession.query(ExperimentLayout).filter(ExperimentLayout.geid == row.geid).first()
            if not layout:
                layout = ExperimentLayout(project=self.project)
                layout.geid = row.geid
                self.dbsession.add(layout)
                self.log.info('Created experiment layout {} in project {}'.format(layout.geid, self.project.geid))
            # Cell line pool
            cell_line = None
            if self.get_value(row.cell_line_name):
                cell_line = self.dbsession.query(CellLine).filter(CellLine.name == row.cell_line_name).first()
                if not cell_line:
                    raise LoaderException('Cell line {} not found'.format(row.cell_line_name))
            # Clone
            if self.get_value(row.clone_name):
                if not cell_line:
                    raise LoaderException('Cannot have a clone without a cell line on row {}'.format(i))
                clone = self.dbsession.query(Clone).filter(Clone.name == str(row.clone_name))\
                                                 .filter(Clone.cell_pool == self.get_string(row.cell_pool))\
                                                 .filter(Clone.project == self.project)\
                                                 .first()
                if not clone:
                    clone = Clone(name=str(row.clone_name), cell_pool=self.get_string(row.cell_pool), project=self.project)
                    clone.cell_line = cell_line
                    self.dbsession.add(clone)
                    self.log.info('Created clone {} in cell line {} pool {}'.format(clone.name, clone.cell_line.name, clone.cell_pool))
            # Well content
            if clone:
                content = WellContent(clone=clone)
                if guide:
                    content.guides.append(guide)
                content.replicate_group = self.get_int(row.replicate_group)
                content.is_control = bool(row.is_control)
                content.content_type = self.to_content_type(row.content_type, i)
                self.dbsession.add(content)
                self.log.info("Created well content for clone {}".format(content.clone.name))
            # Well
            if not row.well_position:
                raise LoaderException('Well position is required on row {}'.format(i))
            well = Well(experiment_layout=layout)
            well.row = row.well_position[0]
            well.column = int(row.well_position[1:])
            well.well_content = content
            self.dbsession.add(well)
            self.log.info('Created well {}{} in layout {}'.format(well.row, well.column, well.experiment_layout.geid))

    # def load_sequencing_libraries(self):
    #     df = self.xls.parse('SequencingLibrary')
    #     for row in df.itertuples():
            sequencing_library = self.dbsession.query(SequencingLibrary).filter(SequencingLibrary.slxid == row.slxid).first()
            if not sequencing_library:
                sequencing_library = SequencingLibrary()
                sequencing_library.slxid = row.slxid
                sequencing_library.library_type = row.library_type
                self.dbsession.add(sequencing_library)
                self.log.info('Created sequening library {}'.format(sequencing_library.slxid))
            experiment_layout = self.dbsession.query(ExperimentLayout).filter(ExperimentLayout.geid == row.experiment_layout_geid).first()
            well = self.dbsession.query(Well).filter(Well.row == row.well_position[0]).filter(Well.column == row.well_position[1:]).filter(Well.experiment_layout_id == experiment_layout.id).first()
            seq_lib_content = SequencingLibraryContent(sequencing_library=sequencing_library, \
                                    well=well, \
                                    dna_source=row.dna_source, \
                                    sequencing_barcode=row.sequencing_barcode, \
                                    sequencing_sample_name = row.sequencing_sample_name)
            self.dbsession.add(seq_lib_content)
            self.log.info('Created sequencing library content {} for library {} in layout {}'.format(seq_lib_content.sequencing_barcode, sequencing_library.slxid, experiment_layout.geid))

    def load_plates(self):
        sheet = self.xls.parse('Plate')
        mandatory_fields = ['layout_geid',
                            'plate_name']
        if not sheet.empty:
            self.check_mandatory_fields('Plate', sheet, mandatory_fields)
        for i, row in enumerate(sheet.itertuples(), 1):
            layout = self.dbsession.query(ExperimentLayout)\
                                   .filter(ExperimentLayout.geid == row.layout_geid).first()
            if not layout:
                raise LoaderException('Layout {} not found (Plate tab, row {})'.format(row.layout_geid, i))
            plate = Plate(experiment_layout=layout)
            plate.name = row.plate_name
            plate.barcode = row.plate_barcode
            plate.description = row.description
            self.dbsession.add(plate)
            self.log.info('Plate {} in layout {} created'.format(plate.plate_name, plate.experiment_layout.geid))


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
                                           .filter(Well.experiment_layout == self.plate.experiment_layout)\
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
                                           .filter(Well.experiment_layout == self.plate.experiment_layout)\
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


# --------------------------------------------------------------------------------
# VariantLoader class
# --------------------------------------------------------------------------------
class VariantLoader(Loader):

    def __init__(self, dbsession, project_geid, workbook_file, variant_caller):
        self.log = logging.getLogger(__name__)
        self.dbsession = dbsession
        self.project_geid = project_geid
        self.xls = pandas.ExcelFile(workbook_file)
        self.variant_caller = variant_caller

    def clean(self):
        variants = self.dbsession.query(VariantResult).\
                   join(SequencingLibraryContent).\
                   join(Well).\
                   join(ExperimentLayout).\
                   join(Project).\
                   filter(Project.geid == self.project_geid).\
                   all()

        for v in variants:
            self.dbsession.delete(v)

        self.dbsession.flush()

        self.log.info('Deleted {} variant results'.format(len(variants)))

    def load_sheet(self, sheet_name, variant_type):
        sheet = self.xls.parse(sheet_name, header=None, skiprows=1)
        header = ['sample',
                  'barcode',
                  'variant_type_consequence',
                  'symbol_gene_id',
                  'cdna_effect',
                  'protein_effect',
                  'codons',
                  'chromosome',
                  'position',
                  'ref',
                  'alt',
                  'igv_links',
                  'specific',
                  'confidence',
                  'allele_fraction',
                  'depth',
                  'filter',
                  'quality',
                  'amplicon',
                  'pubmed',
                  'gene',
                  'exon',
                  'intron',
                  'existing_variation',
                  'sift',
                  'polyphen',
                  'clinical_significance',
                  'gmaf',
                  'offset_from_primer_end',
                  'indel_length',
                  'forward_context',
                  'alleles',
                  'reverse_context',
                  'position_noise_threshold',
                  'library_noise_threshold',
                  'allele_fraction_indep_computed',
                  'depth_indep_computed']
        if len(sheet) == 0:
            return
        sheet.columns = header
        for i, row in enumerate(sheet.itertuples(), 1):
            seq_lib_content = self.dbsession.query(SequencingLibraryContent) \
                                            .filter(SequencingLibraryContent.sequencing_sample_name == row.sample) \
                                            .filter(SequencingLibraryContent.sequencing_barcode == row.barcode) \
                                            .first()
            if not seq_lib_content:
                raise LoaderException("There is no sequencing library content for {} sample with {} barcode".format(row.sample, row.barcode))
            result = VariantResult(sequencing_library_content=seq_lib_content)
            result.variant_caller = self.variant_caller
            result.variant_type = variant_type
            result.consequence = self.get_value(row.variant_type_consequence)
            result.gene_id = self.get_value(row.symbol_gene_id)
            result.cdna_effect = self.get_value(row.cdna_effect)
            result.protein_effect = self.get_value(row.protein_effect)
            result.codons = self.get_value(row.codons)
            result.chromosome = self.get_value(row.chromosome)
            result.position = self.get_int(row.position)
            result.ref = self.get_value(row.ref)
            result.alt = self.get_value(row.alt)
            result.allele_fraction = self.get_float(row.allele_fraction)
            result.depth = self.get_int(row.depth)
            result.quality = self.get_int(row.quality)
            result.amplicon = self.get_value(row.amplicon)
            result.gene = self.get_value(row.gene)
            result.exon = self.get_value(row.exon)
            result.intron = self.get_value(row.intron)
            result.existing_variation = self.get_value(row.existing_variation)
            result.sift = self.get_value(row.sift)
            result.polyphen = self.get_value(row.polyphen)
            result.clinical_significance = self.get_value(row.clinical_significance)
            result.gmaf = self.get_value(row.gmaf)
            result.offset_from_primer_end = self.get_int(row.offset_from_primer_end)
            result.indel_length = self.get_int(row.indel_length)
            result.forward_context = self.get_value(row.forward_context)
            result.alleles = self.get_value(row.alleles)
            result.reverse_context = self.get_value(row.reverse_context)
            self.dbsession.add(result)
            self.log.info("Variant result added for {} sample with {} barcode".format(row.sample, row.barcode))

    def load(self):
        self.load_sheet('SNVs', 'SNV')
        self.load_sheet('Indels', 'INDEL')
