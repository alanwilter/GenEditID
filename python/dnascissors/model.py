"""
Python3 file of the genome-editing project
Created by Anne Pajon @pajanne on 08/03/2017
"""
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship
from sqlalchemy.orm import backref
from sqlalchemy import Table, Column, Integer, String, Date, DateTime, Float, Boolean, Enum, ForeignKey, UniqueConstraint

Base = declarative_base()


guide_well_content_association = Table('guide_well_content_association', Base.metadata,
    Column('well_content_id', Integer, ForeignKey('well_content.id', ondelete='CASCADE')),
    Column('guide_id', Integer, ForeignKey('guide.id', ondelete='CASCADE'))
)


class Project(Base):
    """
    template file: data/templates/YYYYMMDD_GEPXXXXX.xlsx | sheet: Project
    columns: geid	name	scientist	affiliation	group	group_leader	start_date	description
    """
    __tablename__ = 'project'
    id = Column(Integer, primary_key=True)
    geid = Column(String(8), unique=True, nullable=False, index=True)
    name = Column(String(64), unique=True, nullable=False, index=True)
    scientist = Column(String(64))
    group = Column(String(64))
    group_leader = Column(String(64))
    affiliation = Column(String(64))
    start_date = Column(Date, nullable=False)
    end_date = Column(Date, nullable=True)
    project_type = Column(Enum('knock-in', 'knock-out', name='project_type'), nullable=False)
    description = Column(String(1024))
    comments = Column(String(1024))


class Genome(Base):
    __tablename__ = 'genome'
    id = Column(Integer, primary_key=True)
    species = Column(String(32), nullable=False, index=True)
    assembly = Column(String(32), nullable=False, index=True)
    UniqueConstraint('species', 'assembly', name='unique_genome_species_assembly')


class Target(Base):
    """
    template file: data/templates/YYYYMMDD_GEPXXXXX.xlsx | sheet: Target
    columns: project_geid	name	species	assembly	gene_id	chromosome	start	end	strand	description
    """
    __tablename__ = 'target'
    id = Column(Integer, primary_key=True)
    project_id = Column(Integer, ForeignKey('project.id', name='target_project_fk', ondelete='CASCADE'))
    project = relationship(
        Project,
        backref=backref('targets', uselist=True, cascade='delete,all'))
    genome_id = Column(Integer, ForeignKey('genome.id', name='amplicon_genome_fk'), nullable=False)
    genome = relationship(Genome)
    name = Column(String(32), unique=True, nullable=False, index=True)
    gene_id = Column(String(32), nullable=False, index=True)
    chromosome = Column(String(32), nullable=False, index=True)
    start = Column(Integer, nullable=False)
    end = Column(Integer, nullable=False)
    strand = Column(Enum('forward', 'reverse', name='strand'), nullable=False, index=True)
    description = Column(String(1024))


class Guide(Base):
    """
    template file: data/templates/YYYYMMDD_GEPXXXXX.xlsx | sheet: Guide
    columns: target_name	name	guide_sequence	pam_sequence	activity	exon	nuclease
    """
    __tablename__ = 'guide'
    id = Column(Integer, primary_key=True)
    target_id = Column(Integer, ForeignKey('target.id', name='guide_target_fk', ondelete='CASCADE'))
    target = relationship(
        Target,
        backref=backref('guides', uselist=True, cascade='delete,all'))
    genome_id = Column(Integer, ForeignKey('genome.id', name='amplicon_genome_fk'), nullable=False)
    genome = relationship(Genome)
    well_contents = relationship('WellContent', secondary=guide_well_content_association, cascade="all", passive_deletes=True)
    name = Column(String(32), nullable=False, index=True)
    guide_sequence = Column(String(250), nullable=False)
    pam_sequence = Column(String(6), nullable=False)
    activity = Column(Integer, nullable=False)
    exon = Column(Integer, nullable=False)
    nuclease = Column(String(250), nullable=False)


class GuideMismatch(Base):
    """
    template file: data/templates/YYYYMMDD_GEPXXXXX.xlsx | sheet: GuideMismatches
    columns: guide_name	is_off_target_coding_region	number_of_mismatches	number_of_offtargets
    """
    __tablename__ = 'guide_mismatch'
    id = Column(Integer, primary_key=True)
    guide_id = Column(Integer, ForeignKey('guide.id', name='guide_mismatch_guide_fk', ondelete='CASCADE'))
    guide = relationship(
        Guide,
        backref=backref('guide_mismatches', uselist=True, cascade='delete,all'))
    is_off_target_coding_region = Column(Boolean, nullable=False)
    number_of_mismatches = Column(Integer, nullable=False)
    number_of_off_targets = Column(Integer, nullable=False)


class Amplicon(Base):
    __tablename__ = 'amplicon'
    id = Column(Integer, primary_key=True)
    genome_id = Column(Integer, ForeignKey('genome.id', name='amplicon_genome_fk'), nullable=False)
    genome = relationship(Genome)
    dna_feature = Column(Enum('gene', 'precursor', 'non-coding', name='dna_feature'))
    chromosome = Column(String(32), nullable=False, index=True)
    start = Column(Integer, nullable=False)
    end = Column(Integer, nullable=False)


class AmpliconSelection(Base):
    """
    template file: data/templates/YYYYMMDD_GEPXXXXX.xlsx | sheet: AmpliconSelection
    columns: guide_name	experiment_type	guide_location	guide_strand	score	amplicon_name	is_on_target	dna_feature	chromosome	primer_geid	primer_sequence	primer_strand	primer_start	primer_end	description
    """
    __tablename__ = 'amplicon_selection'
    id = Column(Integer, primary_key=True)
    guide_id = Column(Integer, ForeignKey('guide.id', name="ampliconselection_guide_fk", ondelete='CASCADE'))
    guide = relationship(
        Guide,
        backref=backref('amplicon_selections', uselist=True, cascade='delete,all'))
    amplicon_id = Column(Integer, ForeignKey('amplicon.id', name="ampliconselection_amplicon_fk", ondelete='CASCADE'))
    amplicon = relationship(
        Amplicon,
        backref=backref('amplicon_selections', uselist=True, cascade='delete,all'))
    experiment_type = Column(String(32))
    guide_location = Column(Integer, nullable=False)
    guide_strand = Column(Enum('forward', 'reverse', name='strand'), nullable=False)
    is_on_target = Column(Boolean, nullable=False)
    score = Column(Integer)
    description = Column(String(1024))


primer_amplicon_association = Table('primer_amplicon_association', Base.metadata,
    Column('amplicon_id', Integer, ForeignKey('amplicon.id', ondelete='CASCADE')),
    Column('primer_id', Integer, ForeignKey('primer.id', ondelete='CASCADE'))
)


class Primer(Base):
    __tablename__ = 'primer'
    id = Column(Integer, primary_key=True)
    genome_id = Column(Integer, ForeignKey('genome.id', name='amplicon_genome_fk'), nullable=False)
    genome = relationship(Genome)
    amplicons = relationship('Amplicon', secondary=primer_amplicon_association)
    geid = Column(String(8), unique=True, nullable=False, index=True)
    sequence = Column(String(250), nullable=False)
    strand = Column(Enum('forward', 'reverse', name='strand'), nullable=False)
    start = Column(Integer, nullable=False)
    end = Column(Integer, nullable=False)


class CellLine(Base):
    __tablename__ = 'cell_line'
    id = Column(Integer, primary_key=True)
    genome_id = Column(Integer, ForeignKey('genome.id', name='amplicon_genome_fk'), nullable=False)
    genome = relationship(Genome)
    name = Column(String(32), unique=True, nullable=False, index=True)
    pool = Column(String(32))
    description = Column(String(1024))


class Clone(Base):
    __tablename__ = 'clone'
    id = Column(Integer, primary_key=True)
    cell_line_id = Column(Integer, ForeignKey('cell_line.id', name='clone_cellline_fk', ondelete='CASCADE'))
    cell_line = relationship(
        CellLine,
        backref=backref('clones', uselist=True, cascade='delete,all'))
    name = Column(String(32), nullable=False, unique=True, index=True)
    description = Column(String(1024))


class ExperimentLayout(Base):
    """
    template file: data/templates/YYYYMMDD_GEPXXXXX.xlsx | sheet: ExperimentLayout
    columns: project_geid	geid	well_position	cell_line_name	clone_name	guide_name	replicate_group	is_control	content_type
    """
    __tablename__ = 'experiment_layout'
    id = Column(Integer, primary_key=True)
    project_id = Column(Integer, ForeignKey('project.id', name='experiment_layout_project_fk', ondelete='CASCADE'))
    project = relationship(
        Project,
        backref=backref('experiment_layouts', uselist=True, cascade='delete,all'))
    geid = Column(String(12), unique=True, nullable=False, index=True)


class Plate(Base):
    """
    template file: data/templates/YYYYMMDD_GEPXXXXX.xlsx | sheet: Plate
    columns: experiment_layout_geid	plate_barcode	geid	description
    """
    __tablename__ = 'plate'
    id = Column(Integer, primary_key=True)
    experiment_layout_id = Column(Integer, ForeignKey('experiment_layout.id', name='plate_experiment_layout_fk', ondelete='CASCADE'))
    experiment_layout = relationship(
        ExperimentLayout,
        backref=backref('plates', uselist=True, cascade='delete,all'))
    barcode = Column(String(32), index=True)
    geid = Column(String(32), unique=True, nullable=False, index=True)
    description = Column(String)


class WellContent(Base):
    __tablename__ = 'well_content'
    id = Column(Integer, primary_key=True) # is this id the replicate goup?
    clone_id = Column(Integer, ForeignKey('clone.id', name='well_content_clone_fk', ondelete='CASCADE'))
    clone = relationship(
        Clone,
        backref=backref('well_contents', uselist=True, cascade='delete,all'))
    guides = relationship('Guide', secondary=guide_well_content_association, cascade="all", passive_deletes=True)
    replicate_group = Column(Integer, nullable=False, default=0)
    content_type = Column(Enum('wild-type', 'knock-out', 'background', 'normalisation', 'sample', 'empty-vector', 'empty', name='content_type'), nullable=False)
    is_control = Column(Boolean, nullable=False, default=False)


class Well(Base):
    __tablename__ = 'well'
    id = Column(Integer, primary_key=True)
    experiment_layout_id = Column(Integer, ForeignKey('experiment_layout.id', name='well_experiment_layout_fk', ondelete='CASCADE'), nullable=False)
    experiment_layout = relationship(
        ExperimentLayout,
        backref=backref('wells', uselist=True, cascade='delete,all'))
    well_content_id = Column(Integer, ForeignKey('well_content.id', name='well_well_content_fk', ondelete='CASCADE'), nullable=True)
    well_content = relationship(
        WellContent,
        backref=backref('wells', uselist=True, cascade='delete,all'))
    row = Column(String(1), nullable=False)
    column = Column(Integer, nullable=False)
    UniqueConstraint('experiment_layout_id', 'row', 'column', name='well_unique_in_layout')

    def is_empty(self):
        return self.well_content.clone == None


class SequencingLibrary(Base):
    """
    template file: data/templates/YYYYMMDD_GEPXXXXX.xlsx | sheet: SequencingLibrary
    columns: experiment_layout_geid	well_position	dna_source	slxid	library_type	sequencing_barcode	sequencing_sample_name
    """
    __tablename__ = 'sequencing_library'
    id = Column(Integer, primary_key=True)
    slxid = Column(String(12), unique=True, nullable=False)
    library_type = Column(String(32))


class SequencingLibraryContent(Base):
    __tablename__ = 'sequencing_library_content'
    id = Column(Integer, primary_key=True)
    well_id = Column(Integer, ForeignKey('well.id', name='well_sequencing_library_content_fk', ondelete='CASCADE'), nullable=False)
    well = relationship(
        Well,
        backref=backref('sequencing_library_contents', uselist=True, cascade='delete,all'))
    sequencing_library_id = Column(Integer, ForeignKey('sequencing_library.id', name='well_sequencing_library_fk', ondelete='CASCADE'), nullable=False)
    sequencing_library = relationship(
        SequencingLibrary,
        backref=backref('sequencing_library_contents', uselist=True, cascade='delete,all'))
    dna_source = Column(Enum('fixed cells', 'gDNA', 'non-fixed cells', name='dna_source'), nullable=False)
    sequencing_barcode = Column(String(20), nullable=False)
    sequencing_sample_name = Column(String(32), nullable=True)


class VariantResult(Base):
    __tablename__ = 'variant_result'
    id = Column(Integer, primary_key=True)
    sequencing_library_content_id = Column(Integer, ForeignKey('sequencing_library_content.id', name='sequencing_library_content_variant_result_fk', ondelete='CASCADE'), nullable=False)
    sequencing_library_content = relationship(
        SequencingLibraryContent,
        backref=backref('variant_results', uselist=True, cascade='delete,all'))
    variant_caller = Column(Enum('VarDict', 'HaplotypeCaller', 'FreeBayes', name='variant_caller'), nullable=False, index=True)
    variant_type = Column(Enum('INDEL', 'SNV', name='variant_type'), nullable=False, index=True)
    consequence = Column(String(250))
    gene_id = Column(String(32))
    gene = Column(String(32))
    cdna_effect = Column(String(250))
    protein_effect = Column(String(250))
    codons = Column(String(250))
    chromosome = Column(String(32))
    position = Column(Integer)
    ref = Column(String(250))
    alt = Column(String(250))
    allele_fraction = Column(Float)
    depth = Column(Integer)
    quality = Column(Integer)
    amplicon = Column(String(250))
    exon = Column(String(32))
    intron = Column(String(32))
    existing_variation = Column(String(250))
    sift = Column(String(250))
    polyphen = Column(String(250))
    clinical_significance = Column(String(250))
    gmaf = Column(String(250))
    offset_from_primer_end = Column(Integer)
    indel_length = Column(Integer)
    forward_context = Column(String(250))
    alleles = Column(String(250))
    reverse_context = Column(String(250))

    @property
    def frame(self):
        if self.indel_length is None:
            return None
        frame = self.indel_length % 3
        if frame < 0:
            frame += 3
        return frame


class MutationSummary(Base):
    __tablename__ = 'mutation_summary'
    id = Column(Integer, primary_key=True)
    sequencing_library_content_id = Column(Integer, ForeignKey('sequencing_library_content.id', name='sequencing_library_content_mutation_summary_fk', ondelete='CASCADE'), unique=True, nullable=False)
    sequencing_library_content = relationship(
        SequencingLibraryContent,
        backref=backref('mutation_summaries', uselist=True, cascade='delete,all'))
    zygosity = Column(Enum('homo', 'smut', 'dmut', 'iffy', 'warn', name='zygosity'))
    consequence = Column(String(250))
    has_off_target = Column(Boolean)
    has_frameshift = Column(Boolean)
    variant_caller_presence = Column(Enum('VH', 'V-', '-H', 'V?', name='variant_caller_presence'))
    score = Column(Integer)


class ProteinAbundance(Base):
    __tablename__ = 'abundance'
    id = Column(Integer, primary_key=True)
    well_id = Column(Integer, ForeignKey('well.id', name='abundance_well_fk', ondelete='CASCADE'), nullable=False)
    well = relationship(
        Well,
        backref=backref('abundances', uselist=True, cascade='delete,all'))
    plate_id = Column(Integer, ForeignKey('plate.id', name='abundance_plate_fk', ondelete='CASCADE'), nullable=False)
    plate = relationship(
        Plate,
        backref=backref('abundances', uselist=True, cascade='delete,all'))
    intensity_channel_700 = Column(Float)
    intensity_channel_800 = Column(Float)

    @property
    def ratio_800_700(self):
        if self.intensity_channel_800 is None or self.intensity_channel_700 is None:
            return None
        if self.intensity_channel_700 == 0.0:
            return 0.0
        return self.intensity_channel_800 / self.intensity_channel_700


class CellGrowth(Base):
    __tablename__ = 'growth'
    id = Column(Integer, primary_key=True)
    well_id = Column(Integer, ForeignKey('well.id', name='growth_well_fk', ondelete='CASCADE'), nullable=False)
    well = relationship(
        Well,
        backref=backref('growths', uselist=True, cascade='delete,all'))
    plate_id = Column(Integer, ForeignKey('plate.id', name='abundance_plate_fk', ondelete='CASCADE'), nullable=False)
    plate = relationship(
        Plate,
        backref=backref('growths', uselist=True, cascade='delete,all'))
    timestamp = Column(DateTime, nullable=False)
    hours = Column(Integer, nullable=False)
    confluence_percentage = Column(Float, nullable=False)
