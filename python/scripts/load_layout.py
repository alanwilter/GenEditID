import sqlalchemy
import logging
import pandas

from dnascissors.config import cfg
from dnascissors.model import Base
from dnascissors.model import Project
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
from dnascissors.model import VariantResult
from dnascissors.model import ProteinAbundance
from dnascissors.model import CellGrowth

from excelloader import ExcelLoader
import log as logger


class LayoutLoader(ExcelLoader):

    COLUMN_SEPARATOR = ','

    def __init__(self, session, workbook_file):
        self.log = logging.getLogger('dnascissors')
        self.session = session
        self.xls = pandas.ExcelFile(workbook_file)

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
            raise Exception("Well content type not recognised: {:s} (row {:d})".format(str(value), i))
        else:
            return 'empty'

    def load_all(self):
        self.load_projects()
        self.session.commit()
        self.load_targets()
        self.session.commit()
        self.load_guides()
        self.session.commit()
        self.load_amplicon_selection()
        self.session.commit()
        self.load_experiment_layout()
        self.session.commit()
        self.load_plates()
        self.session.commit()
        self.load_sequencing_libraries()
        self.session.commit()

    def load_projects(self):
        sheet = self.xls.parse('Project')
        for i, row in enumerate(sheet.itertuples(), 1):
            if not row.geid:
                raise Exception('Project identifier is required on row {:d}'.format(i))
            project = self.session.query(Project).filter(Project.geid == row.geid).first()
            if project:
                self.log.info("Already have a project {:s} ({:s})".format(project.geid, project.name))
                return
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
            if not row.name:
                raise Exception('Target name is required on row {:d}'.format(i))
            target = self.session.query(Target).filter(Target.name == row.name).first()
            if target:
                self.log.info("Already have a target called {:s}. It's from the project {:s}.".format(target.name, target.project.name))
                return
            target = Target(project=project)
            target.name = row.name
            target.species = row.species
            target.assembly = row.assembly
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
            guide = Guide(target=target)
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
            if not row.amplicon_name:
                raise Exception('Amplicon name is required on row {:d}'.format(i))
            amplicon = self.session.query(Amplicon).filter(Amplicon.name == row.amplicon_name).first()
            if not amplicon:
                amplicon = Amplicon()
                amplicon.name = row.amplicon_name
                amplicon.is_on_target = bool(row.is_on_target)
                amplicon.dna_feature = self.to_dna_feature(row.dna_feature, i)
                amplicon.chromosome = str(row.chromosome)
                #amplicon.start = ??? # should be calculate
                #amplicon.end = ??? # should be calculate
                self.session.add(amplicon)
                self.log.info('Created amplicon {:s}'.format(amplicon.name))
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
                if not row.score and row.score != 'NA':
                    selection.score = int(row.score)
                self.session.add(selection)
                self.log.info('Created amplicon selection of amplicon {:s} at {:d}'.format(amplicon.name, selection.guide_location))
            # Find or create primer
            if not row.primer_geid:
                raise Exception('Primer GEID is required on row {:d}'.format(i))
            primer = self.session.query(Primer).filter(Primer.geid == row.primer_geid).first()
            if not primer:
                primer = Primer()
                primer.geid = row.primer_geid
                primer.sequence = row.primer_sequence
                primer.strand = self.to_strand(row.primer_strand, i)
                primer.start = int(row.primer_start)
                primer.end = int(row.primer_end)
                primer.description = self.get_string(row.description, 1024)
                self.session.add(primer)
                self.log.info('Created primer {:s}'.format(primer.geid))
            primer.amplicons.append(amplicon)
            self.log.info('Linked amplicon {:s} with primer {:s}'.format(amplicon.name, primer.geid))

    def load_experiment_layout(self):
        sheet = self.xls.parse('ExperimentLayout')
        guide = None
        for i, row in enumerate(sheet.itertuples(), 1):
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
            # Experiment layout
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
                cell_line = self.session.query(CellLine).filter(CellLine.name == row.cell_line_name).first()
                if not cell_line:
                    cell_line = CellLine(name=row.cell_line_name)
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


def main():
    import argparse
    import os
    parser = argparse.ArgumentParser()
    parser.add_argument("--layout", dest="file_layout", action="store", help="The file layout e.g. '20170118_GEP00001.xlsx'", required=True)
    parser.add_argument("--clean", dest="clean_db", action="store_true", default=False, help="Clean database before loading?")
    options = parser.parse_args()

    log = logger.get_custom_logger(os.path.join(os.path.dirname(__file__), 'load_layout.log'))

    engine = sqlalchemy.create_engine(cfg['DATABASE_URI'])
    Base.metadata.bind = engine
    DBSession = sqlalchemy.orm.sessionmaker(bind=engine)
    session = DBSession()

    if options.clean_db:

        delete_count = session.query(CellGrowth).delete()
        log.info('Deleted {:d} cell growths'.format(delete_count))

        delete_count = session.query(ProteinAbundance).delete()
        log.info('Deleted {:d} protein abundances'.format(delete_count))

        delete_count = session.query(VariantResult).delete()
        log.info('Deleted {:d} variant results'.format(delete_count))

        delete_count = session.query(Plate).delete()
        log.info('Deleted {:d} plates'.format(delete_count))

        delete_count = session.query(SequencingLibrary).delete()
        log.info('Deleted {:d} sequencing libraries'.format(delete_count))

        delete_count = session.query(SequencingLibraryContent).delete()
        log.info('Deleted {:d} sequencing library contents'.format(delete_count))

        delete_count = session.query(Primer).delete()
        log.info('Deleted {:d} primers'.format(delete_count))

        delete_count = session.query(AmpliconSelection).delete()
        log.info('Deleted {:d} amplicon selections'.format(delete_count))

        delete_count = session.query(Amplicon).delete()
        log.info('Deleted {:d} amplicons'.format(delete_count))

        delete_count = session.query(Guide).delete()
        log.info('Deleted {:d} guides'.format(delete_count))

        delete_count = session.query(CellLine).delete()
        log.info('Deleted {:d} cell lines'.format(delete_count))

        delete_count = session.query(Clone).delete()
        log.info('Deleted {:d} clones'.format(delete_count))

        delete_count = session.query(WellContent).delete()
        log.info('Deleted {:d} well contents'.format(delete_count))

        delete_count = session.query(Well).delete()
        log.info('Deleted {:d} wells'.format(delete_count))

        delete_count = session.query(ExperimentLayout).delete()
        log.info('Deleted {:d} experiment layouts'.format(delete_count))

        delete_count = session.query(Target).delete()
        log.info('Deleted {:d} targets'.format(delete_count))

        delete_count = session.query(Project).delete()
        log.info('Deleted {:d} projects'.format(delete_count))

    loader = LayoutLoader(session, options.file_layout)

    try:
        loader.load_all()
    except Exception as e:
        log.exception(e)
    finally:
        session.close()

if __name__ == '__main__':
    main()
