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


class Genome(Base):
    # Ref from https://www.ncbi.nlm.nih.gov/grc
    __tablename__ = 'genome'
    id = Column(Integer, primary_key=True)
    species = Column(String(32), nullable=False, index=True)
    assembly = Column(String(32), unique=True, nullable=False, index=True)
    UniqueConstraint('species', 'assembly', name='unique_genome_species_assembly')


class CellLine(Base):
    # Ref from http://web.expasy.org/cellosaurus/
    # Data ftp://ftp.expasy.org/databases/cellosaurus/cellosaurus.txt
    __tablename__ = 'cell_line'
    id = Column(Integer, primary_key=True)
    name = Column(String(32), unique=True, nullable=False, index=True)


class Project(Base):
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

    @property
    def is_abundance_data_available(self):
        if self.experiment_layouts:
            for experiment_layout in self.experiment_layouts:
                if experiment_layout.wells:
                    for well in experiment_layout.wells:
                        if len(well.abundances) > 0:
                            return True
        return False

    @property
    def is_growth_data_available(self):
        if self.experiment_layouts:
            for experiment_layout in self.experiment_layouts:
                if experiment_layout.wells:
                    for well in experiment_layout.wells:
                        if len(well.growths) > 0:
                            return True
        return False

    @property
    def is_variant_data_available(self):
        if self.experiment_layouts:
            for experiment_layout in self.experiment_layouts:
                if experiment_layout.wells:
                    for well in experiment_layout.wells:
                        if well.sequencing_library_contents:
                            for sequencing_library_content in well.sequencing_library_contents:
                                if len(sequencing_library_content.variant_results) > 0:
                                    return True
        return False


class Target(Base):
    __tablename__ = 'target'
    id = Column(Integer, primary_key=True)
    project_id = Column(Integer, ForeignKey('project.id', name='target_project_fk', ondelete='CASCADE'))
    project = relationship(
        Project,
        backref=backref('targets', uselist=True, cascade='delete,all'))
    genome_id = Column(Integer, ForeignKey('genome.id', name='target_genome_fk'), nullable=False)
    genome = relationship(Genome)
    name = Column(String(32), nullable=False, index=True)
    gene_id = Column(String(32), nullable=False, index=True)
    chromosome = Column(String(32), nullable=False, index=True)
    start = Column(Integer, nullable=False)
    end = Column(Integer, nullable=False)
    strand = Column(Enum('forward', 'reverse', name='strand'), nullable=False, index=True)
    description = Column(String(1024))
    UniqueConstraint('name', 'project', name='unique_target_in_project')


class Guide(Base):
    __tablename__ = 'guide'
    id = Column(Integer, primary_key=True)
    target_id = Column(Integer, ForeignKey('target.id', name='guide_target_fk', ondelete='CASCADE'))
    target = relationship(
        Target,
        backref=backref('guides', uselist=True, cascade='delete,all'))
    genome_id = Column(Integer, ForeignKey('genome.id', name='guide_genome_fk'), nullable=False)
    genome = relationship(Genome)
    well_contents = relationship('WellContent', secondary=guide_well_content_association, cascade="all", passive_deletes=True)
    name = Column(String(32), nullable=False, index=True)
    guide_sequence = Column(String(250), nullable=False)
    pam_sequence = Column(String(6), nullable=False)
    activity = Column(Integer, nullable=False)
    exon = Column(Integer, nullable=False)
    nuclease = Column(String(250), nullable=False)


class GuideMismatch(Base):
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
    genome_id = Column(Integer, ForeignKey('genome.id', name='primer_genome_fk'), nullable=False)
    genome = relationship(Genome)
    amplicons = relationship('Amplicon', secondary=primer_amplicon_association)
    geid = Column(String(32), unique=True, nullable=False, index=True)
    sequence = Column(String(250), nullable=False)
    strand = Column(Enum('forward', 'reverse', name='strand'), nullable=False)
    start = Column(Integer, nullable=False)
    end = Column(Integer, nullable=False)


class Clone(Base):
    __tablename__ = 'clone'
    id = Column(Integer, primary_key=True)
    name = Column(String(32), nullable=False, index=True)
    cell_line_id = Column(Integer, ForeignKey('cell_line.id', name='clone_cell_line_fk'))
    cell_line = relationship(CellLine)
    cell_pool = Column(String(32))
    project_id = Column(Integer, ForeignKey('project.id', name='clone_project_fk', ondelete='CASCADE'))
    project = relationship(
        Project,
        backref=backref('clones', uselist=True, cascade='delete,all'))
    description = Column(String(1024))
    UniqueConstraint('name', 'cell_pool', 'project', name='unique_clone_in_project')


class ExperimentLayout(Base):
    __tablename__ = 'experiment_layout'
    id = Column(Integer, primary_key=True)
    project_id = Column(Integer, ForeignKey('project.id', name='experiment_layout_project_fk', ondelete='CASCADE'))
    project = relationship(
        Project,
        backref=backref('experiment_layouts', uselist=True, cascade='delete,all'))
    geid = Column(String(12), unique=True, nullable=False, index=True)


class Plate(Base):
    __tablename__ = 'plate'
    id = Column(Integer, primary_key=True)
    experiment_layout_id = Column(Integer, ForeignKey('experiment_layout.id', name='plate_experiment_layout_fk', ondelete='CASCADE'))
    experiment_layout = relationship(
        ExperimentLayout,
        backref=backref('plates', uselist=True, cascade='delete,all'))
    barcode = Column(String(32), index=True)
    geid = Column(String(32), unique=True, nullable=False, index=True)
    description = Column(String)

    @property
    def is_abundance_plate(self):
        if len(self.abundances) > 0:
            return True

    @property
    def is_growth_plate(self):
        if len(self.growths) > 0:
            return True

    @property
    def is_ngs_plate(self):
        for well in self.experiment_layout.wells:
            if well.sequencing_library_contents:
                for sequencing_library_content in well.sequencing_library_contents:
                    if len(sequencing_library_content.variant_results) > 0:
                        return True

    @property
    def plate_type(self):
        types = {'abundance': self.is_abundance_plate,
                 'growth': self.is_growth_plate,
                 #'NGS': self.is_ngs_plate
                 }
        return ','.join([k for k, v in types.items() if v is True])


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
        if self.indel_length:
            if abs(self.indel_length) > 3:
                return abs(self.indel_length) % 3
            return abs(self.indel_length)

    @property
    def variant_str_summary(self):
        return "{}:{}:{}:{}:{:.3f}:{}".format(self.variant_caller,
                                              self.variant_type,
                                              self.frame,
                                              self.gene_id,
                                              self.allele_fraction,
                                              self.alleles)


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
