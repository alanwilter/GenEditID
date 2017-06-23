import logging
from datetime import datetime
import csv
import pandas

from dnascissors.model import Project
from dnascissors.model import Genome
from dnascissors.model import Target
from dnascissors.model import Guide
from dnascissors.model import Amplicon
from dnascissors.model import AmpliconSelection
from dnascissors.model import Primer
from dnascissors.model import CellLine
from dnascissors.model import Clone
from dnascissors.model import ExperimentLayout
from dnascissors.model import Plate
from dnascissors.model import WellContent
from dnascissors.model import Well
from dnascissors.model import SequencingLibrary
from dnascissors.model import SequencingLibraryContent
from dnascissors.model import CellGrowth
from dnascissors.model import ProteinAbundance
from dnascissors.model import VariantResult
from dnascissors.model import MutationSummary


# --------------------------------------------------------------------------------
# Loader exceptions
# --------------------------------------------------------------------------------
class ExistingEntityException(Exception):

    def __init__(self, objtype, key, msg=None):
        if msg is None:
            msg="Already have a {:s} identified by {:s}.".format(type(objtype).__name__, str(key))

        super(ExistingEntityException, self).__init__(msg)

        self.object_type = type
        self.key = key
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
            raise ValueError('Strand must be "forward" or "reverse", not "{:s}" (row {:d})'.format(str(value), i))

    def to_dna_feature(self, value, i):
        if self.get_value(value):
            if value in ['gene', 'precursor', 'non-coding']:
                return value
            raise ValueError('DNA feature must be "gene", "precursor", "non-coding", not "{:s}" (row {:d})'.format(str(value), i))

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
            raise Exception("Well content type not recognised: {:s} (row {:d})".format(str(value), i))
        else:
            return 'empty'


# --------------------------------------------------------------------------------
# LayoutLoader class
# --------------------------------------------------------------------------------
class LayoutLoader(Loader):

    def __init__(self, session, workbook_file):
        self.log = logging.getLogger('dnascissors')
        self.session = session
        self.xls = pandas.ExcelFile(workbook_file)
        self.genome = None

    def clean(self, project):
        self.session.delete(project)
        self.session.flush()

    def load_all(self, clean_if_exists=False):
        self.load_projects(clean_if_exists)
        self.load_targets()
        self.load_guides()
        self.load_amplicon_selection()
        self.load_experiment_layout()
        self.load_plates()
        self.load_sequencing_libraries()

    def load_projects(self, clean_if_exists=False):
        sheet = self.xls.parse('Project')
        for i, row in enumerate(sheet.itertuples(), 1):
            if not row.geid:
                raise Exception('Project identifier is required on row {:d}'.format(i))
            project = self.session.query(Project).filter(Project.geid == row.geid).first()
            if project:
                if clean_if_exists:
                    self.log.info("Already have a project {:s} ({:s})".format(project.geid, project.name))
                    self.log.info("Removing this project and its associated data.")
                    self.clean(project)
                else:
                    raise ExistingEntityException(Project, row.geid, "Already have project {:s}. Will not overwrite it.".format(project.geid))

            project = Project()
            project.geid = row.geid
            project.name = row.name
            project.scientist = row.scientist
            project.institute = row.institute
            project.group = row.group
            project.group_leader = row.group_leader
            project.start_date = self.get_date(row.start_date)
            project.description = self.get_string(row.description, 1024)
            self.session.add(project)
            self.log.info('Created project {:s}'.format(project.name))

    def load_targets(self):
        sheet = self.xls.parse('Target')
        for i, row in enumerate(sheet.itertuples(), 1):
            if not row.project_geid:
                raise Exception('Project identifier is required on row {:d}'.format(i))
            project = self.session.query(Project).filter(Project.geid == row.project_geid).first()
            if not project:
                raise Exception('Project "{:s}" does not exist (row {:d})'.format(row.project_geid, i))
            if not row.species:
                raise Exception('Species is required on row {:d}'.format(i))
            if not row.assembly:
                raise Exception('Genome assembly is required on row {:d}'.format(i))
            # Find or create genome
            genome = self.session.query(Genome).filter(Genome.species == row.species).filter(Genome.assembly == row.assembly).first()
            if not genome:
                genome = Genome(species=row.species, assembly=row.assembly)
                self.session.add(genome)
                self.log.info('Created genome species {:s}, assembly {:s}'.format(genome.species, genome.assembly))
            self.genome = genome
            # Find or create target
            if not row.name:
                raise Exception('Target name is required on row {:d}'.format(i))
            target = self.session.query(Target).filter(Target.name == row.name).first()
            if target:
                self.log.info("Already have a target called {:s}. It's from the project {:s}.".format(target.name, target.project.name))
                return
            target = Target(project=project, genome=self.genome)
            target.name = row.name
            target.gene_id = row.gene_id
            target.chromosome = str(row.chromosome)
            target.start = int(row.start)
            target.end = int(row.end)
            target.strand = self.to_strand(row.strand, i)
            target.description = self.get_string(row.description, 1024)
            self.session.add(target)
            self.log.info('Created target {:s}'.format(target.name))

    def load_guides(self):
        sheet = self.xls.parse('Guide')
        for i, row in enumerate(sheet.itertuples(), 1):
            if not row.target_name:
                raise Exception('Target name is required on row {:d}'.format(i))
            target = self.session.query(Target).filter(Target.name == row.target_name).first()
            if not target:
                raise Exception('Target "{:s}" does not exist (row {:d})'.format(row.target_name, i))
            if not row.name:
                raise Exception('Guide name is required on row {:d}'.format(i))
            guide = Guide(target=target, genome=self.genome)
            guide.name = row.name
            guide.guide_sequence = row.guide_sequence
            guide.pam_sequence = row.pam_sequence
            guide.activity = int(row.activity)
            guide.exon = int(row.exon)
            guide.nuclease = row.nuclease
            self.session.add(guide)
            self.log.info('Created guide {:s}'.format(guide.name))

    def load_amplicon_selection(self):
        sheet = self.xls.parse('AmpliconSelection')
        previous_strand = None
        for i, row in enumerate(sheet.itertuples(), 1):
            # Need the guide first.
            if not row.guide_name:
                raise Exception('Guide name is required on row {:d}'.format(i))
            guide = self.session.query(Guide).filter(Guide.name == row.guide_name).first()
            if not guide:
                raise Exception('Guide "{:s}" does not exist (row {:d})'.format(row.guide_name, i))
            # Find or create the amplicon
            amplicon = self.session.query(Amplicon).filter(Amplicon.genome == self.genome)\
                                                   .filter(Amplicon.chromosome == str(row.chromosome))\
                                                   .filter(Amplicon.start == int(row.forward_primer_start))\
                                                   .filter(Amplicon.end == int(row.reverse_primer_end)).first()
            if not amplicon:
                amplicon = Amplicon(genome=self.genome)
                amplicon.dna_feature = self.to_dna_feature(row.dna_feature, i)
                amplicon.chromosome = str(row.chromosome)
                amplicon.start = int(row.forward_primer_start)
                amplicon.end = int(row.reverse_primer_end)
                self.session.add(amplicon)
                self.log.info('Created amplicon {:s}_chr{:s}_{:d}'.format(amplicon.genome.assembly, amplicon.chromosome, amplicon.start))
            # Find or create the amplicon selection
            strand = self.to_strand(row.guide_strand, i)
            if strand:
                previous_strand = strand
            else:
                if previous_strand:
                    strand = previous_strand
                else:
                    raise Exception("Have no guide strand on row {:d}".format(i))
            selection = self.session.query(AmpliconSelection)\
                                    .filter(AmpliconSelection.amplicon == amplicon)\
                                    .filter(AmpliconSelection.guide_location == int(row.guide_location))\
                                    .filter(AmpliconSelection.guide_strand == strand).first()
            if not selection:
                selection = AmpliconSelection(guide=guide, amplicon=amplicon)
                selection.experiment_type = row.experiment_type
                selection.guide_location = int(row.guide_location)
                selection.guide_strand = strand
                selection.is_on_target = bool(row.is_on_target)
                if pandas.notnull(row.score): # and row.score != 'NA':
                    selection.score = int(row.score)
                selection.description = self.get_string(row.description, 1024)
                self.session.add(selection)
                self.log.info('Created amplicon selection of amplicon {:s}_chr{:s}_{:d} at {:d}'.format(amplicon.genome.assembly, amplicon.chromosome, amplicon.start, selection.guide_location))
            # Find or create forward primer
            if not row.forward_primer_geid:
                raise Exception('Forward primer GEID is required on row {:d}'.format(i))
            forward_primer = self.session.query(Primer).filter(Primer.geid == row.forward_primer_geid).first()
            if not forward_primer:
                forward_primer = Primer(genome=self.genome)
                forward_primer.geid = row.forward_primer_geid
                forward_primer.sequence = row.forward_primer_sequence
                forward_primer.strand = 'forward'
                forward_primer.start = int(row.forward_primer_start)
                forward_primer.end = int(row.forward_primer_end)
                self.session.add(forward_primer)
                self.log.info('Created primer {:s}'.format(forward_primer.geid))
            forward_primer.amplicons.append(amplicon)
            self.log.info('Linked amplicon {:s}_chr{:s}_{:d} with forward primer {:s}'.format(amplicon.genome.assembly, amplicon.chromosome, amplicon.start, forward_primer.geid))
            # Find or create reverse primer
            if not row.reverse_primer_geid:
                raise Exception('Reverse primer GEID is required on row {:d}'.format(i))
            reverse_primer = self.session.query(Primer).filter(Primer.geid == row.reverse_primer_geid).first()
            if not reverse_primer:
                reverse_primer = Primer(genome=self.genome)
                reverse_primer.geid = row.reverse_primer_geid
                reverse_primer.sequence = row.reverse_primer_sequence
                reverse_primer.strand = 'reverse'
                reverse_primer.start = int(row.reverse_primer_start)
                reverse_primer.end = int(row.reverse_primer_end)
                self.session.add(reverse_primer)
                self.log.info('Created primer {:s}'.format(reverse_primer.geid))
            reverse_primer.amplicons.append(amplicon)
            self.log.info('Linked amplicon {:s}_chr{:s}_{:d} with forward primer {:s}'.format(amplicon.genome.assembly, amplicon.chromosome, amplicon.start, reverse_primer.geid))

    def load_experiment_layout(self):
        sheet = self.xls.parse('ExperimentLayout')
        for i, row in enumerate(sheet.itertuples(), 1):
            guide = None
            clone = None
            content = None
            # Need to find project.
            if not row.project_geid:
                raise Exception('Project identifier is required on row {:d}'.format(i))
            project = self.session.query(Project).filter(Project.geid == row.project_geid).first()
            if not project:
                raise Exception('Project "{:s}" does not exist (row {:d})'.format(row.project_geid, i))
            # May also need to find a guide.
            if self.get_value(row.guide_name):
                guide = self.session.query(Guide).filter(Guide.name == row.guide_name).first()
                if not guide:
                    raise Exception("There is no guide with the name {:s}".format(str(row.guide_name)))
            # Find or create experiment layout
            if not row.geid:
                raise Exception('Experiment layout GEID is required on row {:d}'.format(i))
            layout = self.session.query(ExperimentLayout).filter(ExperimentLayout.geid == row.geid).first()
            if not layout:
                layout = ExperimentLayout(project=project)
                layout.geid = row.geid
                self.session.add(layout)
                self.log.info('Created experiment layout {:s} in project {:s}'.format(layout.geid, project.geid))
            # Cell line
            if self.get_value(row.cell_line_name):
                cell_line = self.session.query(CellLine).filter(CellLine.name == row.cell_line_name).filter(CellLine.genome == self.genome).first()
                if not cell_line:
                    cell_line = CellLine(name=row.cell_line_name, genome=self.genome)
                    self.session.add(cell_line)
                    self.log.info('Created cell line {:s}'.format(cell_line.name))
            # Clone
            if self.get_value(row.clone_name):
                if not cell_line:
                    raise Exception('Cannot have a clone without a cell line on row {:d}'.format(i))
                clone = self.session.query(Clone).filter(Clone.cell_line == cell_line).filter(Clone.name == row.clone_name).first()
                if not clone:
                    clone = Clone(cell_line=cell_line)
                    clone.name = row.clone_name
                    self.session.add(clone)
                    self.log.info('Created clone {:s} from cell row {:s}'.format(clone.name, clone.cell_line.name))
            # Well content
            if clone:
                content = WellContent(clone=clone)
                if guide:
                    content.guides.append(guide)
                content.replicate_group = self.get_int(row.replicate_group)
                content.is_control = bool(row.is_control)
                content.content_type = self.to_content_type(row.content_type, i)
                self.session.add(content)
                self.log.info("Created well content for clone {:s}".format(content.clone.name))
            # Well
            if not row.well_position:
                raise Exception('Well position is required on row {:d}'.format(i))
            well = Well(experiment_layout=layout)
            well.row = row.well_position[0]
            well.column = int(row.well_position[1:])
            well.well_content = content
            self.session.add(well)
            self.log.info('Created well {:s}{:d} in layout {:s}'.format(well.row, well.column, well.experiment_layout.geid))

    def load_plates(self):
        df = self.xls.parse('Plate')
        for row in df.itertuples():
            experiment_layout = self.session.query(ExperimentLayout).filter(ExperimentLayout.geid == row.experiment_layout_geid).first()
            if not experiment_layout:
                raise Exception('Experiment layout GEID {:s} is required for plate GEID {:s}'.format(row.experiment_layout_geid, row.geid))
            plate = Plate(experiment_layout=experiment_layout)
            plate.geid = row.geid
            plate.barcode = row.plate_barcode
            plate.description = row.description
            self.session.add(plate)
            self.log.info('Created plate {:s} in layout {:s}'.format(plate.geid, plate.experiment_layout.geid))

    def load_sequencing_libraries(self):
        df = self.xls.parse('SequencingLibrary')
        for row in df.itertuples():
            sequencing_library = self.session.query(SequencingLibrary).filter(SequencingLibrary.slxid == row.slxid).first()
            if not sequencing_library:
                sequencing_library = SequencingLibrary()
                sequencing_library.slxid = row.slxid
                sequencing_library.library_type = row.library_type
                sequencing_library.barcode_size = int(row.barcode_size)
                self.session.add(sequencing_library)
                self.log.info('Created sequening library {:s}'.format(sequencing_library.slxid))
            experiment_layout = self.session.query(ExperimentLayout).filter(ExperimentLayout.geid == row.experiment_layout_geid).first()
            well = self.session.query(Well).filter(Well.row == row.well_position[0]).filter(Well.column == row.well_position[1:]).filter(Well.experiment_layout_id == experiment_layout.id).first()
            seq_lib_content = SequencingLibraryContent(sequencing_library=sequencing_library, \
                                    well=well, \
                                    dna_source=row.dna_source, \
                                    sequencing_barcode=row.sequencing_barcode, \
                                    sequencing_sample_name = row.sequencing_sample_name)
            self.session.add(seq_lib_content)
            self.log.info('Created sequencing library content {:s} for library {:s} in layout {:s}'.format(seq_lib_content.sequencing_barcode, sequencing_library.slxid, experiment_layout.geid))


# --------------------------------------------------------------------------------
# ProteinAbundanceLoader class
# --------------------------------------------------------------------------------
class ProteinAbundanceLoader(Loader):

    def __init__(self, session, csv_file, plate_id):
        self.log = logging.getLogger('dnascissors')
        self.session = session
        self.csv_file = csv_file
        self.plate_id = plate_id
        self.plate = self.session.query(Plate).filter(Plate.geid == plate_id).first()
        if not self.plate:
            raise Exception('No plate {:s}'.format(plate_id))

    def load(self):
        self.log.info("Loading InCell Western signal for plate {:s} from {:s}".format(self.plate_id, self.csv_file))
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
                        well = self.session.query(Well)\
                                           .filter(Well.experiment_layout == self.plate.experiment_layout)\
                                           .filter(Well.row == rowchar)\
                                           .filter(Well.column == column)\
                                           .first()
                        if not well:
                            raise Exception("There is no well at {:s}{:d} on plate {:s}".format(rowchar, column, self.plate.geid))
                        self.log.info("Setting well {:s}{:d} on {:s} for ICW {:d}".format(well.row, well.column, self.plate.geid, channel))

                        abundance = self.session.query(ProteinAbundance)\
                                                .filter(ProteinAbundance.well == well)\
                                                .filter(ProteinAbundance.plate == self.plate)\
                                                .first()
                        if not abundance:
                            abundance = ProteinAbundance(well=well, plate=self.plate)
                            self.session.add(abundance)
                        if channel == 700:
                            abundance.intensity_channel_700 = self.get_float(line[column - 1])
                        elif channel == 800:
                            abundance.intensity_channel_800 = self.get_float(line[column - 1])
                    row = row + 1


# --------------------------------------------------------------------------------
# CellGrowthLoader class
# --------------------------------------------------------------------------------
class CellGrowthLoader(Loader):

    def __init__(self, session, csv_file, plate_id):
        self.log = logging.getLogger('dnascissors')
        self.session = session
        self.csv_file = csv_file
        self.plate_id = plate_id
        self.plate = self.session.query(Plate).filter(Plate.geid == plate_id).first()
        if not self.plate:
            raise Exception('No plate {:s}'.format(plate_id))

    def load(self):
        self.log.info("Loading Incucyte growth information for plate {:s} from {:s}".format(self.plate_id, self.csv_file))

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
                    hour = int(line[1])
                    for column in range(2, len(line)):
                        location = header[column]
                        wellrow = location[0]
                        wellcolumn = int(location[1:])
                        well = self.session.query(Well)\
                                           .filter(Well.experiment_layout == self.plate.experiment_layout)\
                                           .filter(Well.row == wellrow)\
                                           .filter(Well.column == wellcolumn)\
                                           .first()
                        if not well:
                            raise Exception("There is no well at {:s}{:d} on plate {:s}".format(wellrow, wellcolumn, self.plate.geid))
                        value = float(line[column])
                        self.log.info("Setting well {:s}{:d} on {:s} for Incucyte {:f} at {:s}".format(well.row, well.column, self.plate.geid, value, line[0]))
                        growth = self.session.query(CellGrowth)\
                                             .filter(CellGrowth.plate == self.plate)\
                                             .filter(CellGrowth.well == well)\
                                             .filter(CellGrowth.hours == hour)\
                                             .first()
                        if not growth:
                            growth = CellGrowth(well=well, plate=self.plate)
                            growth.hours = hour
                            growth.timestamp = time
                            self.session.add(growth)
                        growth.confluence_percentage = value


# --------------------------------------------------------------------------------
# VariantLoader class
# --------------------------------------------------------------------------------
class VariantLoader(Loader):

    def __init__(self, session, project_geid, workbook_file, variant_caller):
        self.log = logging.getLogger('dnascissors')
        self.session = session
        self.project_geid = project_geid
        self.xls = pandas.ExcelFile(workbook_file)
        self.variant_caller = variant_caller

    def clean(self):
        delete_count = self.session.query(VariantResult).delete()
        self.log.info('Deleted {:d} variant results'.format(delete_count))

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
                  'reverse_context']
        sheet.columns = header
        for i, row in enumerate(sheet.itertuples(), 1):
            seq_lib_content = self.session.query(SequencingLibraryContent) \
                                          .filter(SequencingLibraryContent.sequencing_sample_name == row.sample) \
                                          .filter(SequencingLibraryContent.sequencing_barcode == row.barcode) \
                                          .first()
            if not seq_lib_content:
                raise Exception("There is no sequencing library content for {:s} sample with {:s} barcode".format(row.sample, row.barcode))
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
            self.session.add(result)
            self.log.info("Variant result added for {:s} sample with {:s} barcode".format(row.sample, row.barcode))

    def load(self):
        self.load_sheet('SNVs', 'SNV')
        self.load_sheet('Indels', 'INDEL')


# --------------------------------------------------------------------------------
# MutationLoader class
# --------------------------------------------------------------------------------
class MutationLoader(Loader):

    def __init__(self, session, project_geid):
        self.log = logging.getLogger('dnascissors')
        self.session = session
        self.project_geid = project_geid

    def clean(self):
        delete_count = self.session.query(MutationSummary).delete()
        self.log.info('Deleted {:d} mutation summaries'.format(delete_count))

    def __get_mutations__(self, well, variant_results, variant_caller='VarDict'):
        mutations = []
        nb_variants = 0
        for variant in variant_results:
            # select INDEL variants only from VarDict caller
            if variant.variant_caller == variant_caller and variant.variant_type == 'INDEL':
                nb_variants += 1
                self.log.info('--- variant: {} {} {}'.format(variant.variant_caller, variant.variant_type, variant.allele_fraction))
                # check amplicon and allele are on the same chromosome
                if not variant.amplicon.split('_')[1] == variant.chromosome:
                    raise ValueError('Amplicon on {:s} and variant on {:s}'.format(variant.amplicon.split('_')[1], variant.chromosome))
                if len(well.well_content.guides) == 1:
                    # find which cut sites match the variant
                    guide = well.well_content.guides[0]
                    matched_amplicon_selection = None
                    for amplicon_selection in guide.amplicon_selections:
                        if variant.chromosome.endswith(amplicon_selection.amplicon.chromosome) and variant.position >= amplicon_selection.amplicon.start and variant.position <= amplicon_selection.amplicon.end:
                            matched_amplicon_selection = amplicon_selection
                    if matched_amplicon_selection:
                        # select mutation only if within the amplicon range
                        mutations.append(variant)
                elif len(well.well_content.guides) > 1:
                    raise Exception('More than one associated guide, {} found. Cannot calculate the score.'.format(len(well.well_content.guides)))
            return mutations, nb_variants

    def __characterise_mutations__(self, mutations, nb_variants):
        # caracterise the mutations
        mutation_zygosity = None
        mutation_consequence = None
        mutation_has_off_target = None
        if len(mutations) == nb_variants:
            if len(mutations) > 0:
                mutation_has_off_target = False
                if len(mutations) == 1:
                    if mutations[0].allele_fraction > 0.85:
                        mutation_zygosity = 'homo'
                    elif mutations[0].allele_fraction < 0.85 and mutations[0].allele_fraction > 0.35:
                        mutation_zygosity = 'smut'
                    else:
                        mutation_zygosity = 'iffy'
                elif len(mutations) == 2:
                    if mutations[0].allele_fraction < 0.85 and mutations[0].allele_fraction > 0.35 and mutations[1].allele_fraction < 0.85 and mutations[1].allele_fraction > 0.35:
                        mutation_zygosity = 'dmut'
                    else:
                        mutation_zygosity = 'iffy'
                else:
                    mutation_zygosity = 'iffy'
                mutation_consequence = ','.join(set([m.consequence for m in mutations]))
        else:
            mutation_has_off_target = True
        return mutation_zygosity, mutation_consequence, mutation_has_off_target

    def load(self):
        results = self.session.query(SequencingLibraryContent)\
                              .join(SequencingLibraryContent.variant_results)\
                              .join(SequencingLibraryContent.well)\
                              .join(Well.well_content)\
                              .join(Well.experiment_layout)\
                              .join(ExperimentLayout.project)\
                              .filter(Project.geid == self.project_geid)\
                              .all()
        for sequencing_library_content in results:
            well = sequencing_library_content.well
            layout = well.experiment_layout
            mutations, nb_variants = self.__get_mutations__(well, sequencing_library_content.variant_results, 'VarDict')
            mutation_zygosity, mutation_consequence, mutation_has_off_target = self.__characterise_mutations__(mutations, nb_variants)
            summary = MutationSummary(sequencing_library_content=sequencing_library_content)
            summary.zygosity = mutation_zygosity
            summary.consequence = mutation_consequence
            summary.has_off_target = mutation_has_off_target
            self.session.add(summary)
            self.log.info('[{} {} {}{}] sample: {} | Mutation added: {}\t{}\t{}'.format(layout.project.geid, layout.geid, well.row, well.column, sequencing_library_content.sequencing_barcode, summary.zygosity, summary.has_off_target, summary.consequence))
