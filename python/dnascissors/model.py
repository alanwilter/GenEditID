"""
Python3 file of the genome-editing project
Created by Anne Pajon @pajanne on 08/03/2017
"""
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship
from sqlalchemy.orm import backref
from sqlalchemy import Table, Column, Integer, String, Date, DateTime, Float, Boolean, Enum, ForeignKey, UniqueConstraint

Base = declarative_base()


class Project(Base):
    """
    template file: data/templates/YYYYMMDD_GEPXXXXX_Project.csv
    geid	name	scientist	institute	group	group_leader	start_date	description
    """
    __tablename__ = 'project'
    id = Column(Integer, primary_key=True)
    geid = Column(String(8), unique=True, nullable=False, index=True)
    name = Column(String(64), unique=True, nullable=False, index=True)
    scientist = Column(String(64))
    institute = Column(String(64))
    group = Column(String(64))
    group_leader = Column(String(64))
    start_date = Column(Date, nullable=False)
    end_date = Column(Date, nullable=True)
    description = Column(String(1024))


class Target(Base):
    """
    template file: data/templates/YYYYMMDD_GEPXXXXX_Target.csv
    project_geid	name	species	assembly	gene_id	chromosome	start	end	strand	description
    """
    __tablename__ = 'target'
    id = Column(Integer, primary_key=True)
    project_id = Column(Integer, ForeignKey('project.id', name='target_project_fk', ondelete='cascade'))
    project = relationship(
        Project,
        backref=backref('targets', uselist=True, cascade='delete,all'))
    name = Column(String(32), nullable=False, index=True)
    species = Column(String(32), nullable=False, index=True)
    assembly = Column(String(32), nullable=False, index=True)
    gene_id = Column(String(32), nullable=False, index=True)
    chromosome = Column(String(32), nullable=False, index=True)
    start = Column(Integer, nullable=False)
    end = Column(Integer, nullable=False)
    strand = Column(Enum('forward', 'reverse', name='strand'), nullable=False, index=True)
    description = Column(String(1024))


class Guide(Base):
    """
    template file: data/template/YYYYMMDD_GEPXXXXX_Guide.csv
    target_name	name	guide_sequence	pam_sequence	activity	exon	nuclease
    """
    __tablename__ = 'guide'
    id = Column(Integer, primary_key=True)
    target_id = Column(Integer, ForeignKey('target.id', name='guide_target_fk', ondelete='cascade'))
    target = relationship(
        Target,
        backref=backref('guides', uselist=True, cascade='delete,all'))
    #amplicon_selections = relationship('Amplicon', back_populates="guide")
    name = Column(String(32), nullable=False, index=True)
    guide_sequence = Column(String(250), nullable=False)
    pam_sequence = Column(String(6), nullable=False)
    activity = Column(Integer, nullable=False)
    exon = Column(Integer, nullable=False)
    nuclease = Column(String(250), nullable=False)


class Amplicon(Base):
    __tablename__ = 'amplicon'
    id = Column(Integer, primary_key=True)
    #amplicon_selections = relationship('AmpliconSelection', back_populates="amplicon")
    name = Column(String(32), unique=True, nullable=False, index=True)
    is_on_target = Column(Boolean, nullable=False)
    dna_feature = Column(Enum('gene', 'precursor', 'non-coding', name='dna_feature'))
    chromosome = Column(String(32), nullable=False, index=True)
    start = Column(Integer, nullable=True) # should be calculate when loading primers
    end = Column(Integer, nullable=True) # should be calculate when loading primers


class AmpliconSelection(Base):
    """
    template file: data/templates/YYYYMMDD_GEPXXXXX_AmpliconSelection.csv
    guide_name	experiment_type	guide_location	guide_strand	score	amplicon_name	is_on_target	dna_feature	chromosome	primer_geid	primer_sequence	primer_strand	primer_start	primer_end
    """
    __tablename__ = 'amplicon_selection'
    id = Column(Integer, primary_key=True)
    guide_id = Column(Integer, ForeignKey('guide.id', name="ampliconselection_guide_fk", ondelete='cascade'))
    guide = relationship(
        Guide,
        backref=backref('amplicon_selections', uselist=True, cascade='delete,all'))
    amplicon_id = Column(Integer, ForeignKey('amplicon.id', name="ampliconselection_amplicon_fk", ondelete='cascade'))
    amplicon = relationship(
        Amplicon,
        backref=backref('amplicon_selections', uselist=True, cascade='delete,all'))
    experiment_type = Column(String(32))
    guide_location = Column(Integer, nullable=False)
    guide_strand = Column(Enum('forward', 'reverse', name='strand'), nullable=False)
    score = Column(Integer)


primer_amplicon_association = Table('primer_amplicon_association', Base.metadata,
    Column('amplicon_id', Integer, ForeignKey('amplicon.id', ondelete='cascade')),
    Column('primer_id', Integer, ForeignKey('primer.id', ondelete='cascade'))
)


class Primer(Base):
    __tablename__ = 'primer'
    id = Column(Integer, primary_key=True)
    amplicons = relationship('Amplicon', secondary=primer_amplicon_association)
    geid = Column(String(8), unique=True, nullable=False, index=True)
    sequence = Column(String(250), nullable=False)
    strand = Column(Enum('forward', 'reverse', name='strand'), nullable=False)
    start = Column(Integer, nullable=False)
    end = Column(Integer, nullable=False)
    description = Column(String(1024))


class CellLine(Base):
    __tablename__ = 'cell_line'
    id = Column(Integer, primary_key=True)
    name = Column(String(32), unique=True, nullable=False, index=True)
    description = Column(String(1024))


class Clone(Base):
    __tablename__ = 'clone'
    id = Column(Integer, primary_key=True)
    cell_line_id = Column(Integer, ForeignKey('cell_line.id', name='clone_cellline_fk', ondelete='cascade'))
    cell_line = relationship(
        CellLine,
        backref=backref('clones', uselist=True, cascade='delete,all'))
    name = Column(String(32), nullable=False, unique=True, index=True)
    description = Column(String(1024))


class ExperimentLayout(Base):
    """
    template file: data/templates/YYYYMMDD_GEPXXXXX_ExperimentLayout.csv
    geid	well_position	cell_line_name	clone_name	guide_name	replicate_group	content_type	is_control
    """
    __tablename__ = 'experiment_layout'
    id = Column(Integer, primary_key=True)
    project_id = Column(Integer, ForeignKey('project.id', name='experiment_layout_project_fk', ondelete='cascade'))
    project = relationship(
        Project,
        backref=backref('experiment_layouts', uselist=True, cascade='delete,all'))
    geid = Column(String(12), unique=True, nullable=False, index=True)


class Plate(Base):
    __tablename__ = 'plate'
    id = Column(Integer, primary_key=True)
    experiment_layout_id = Column(Integer, ForeignKey('experiment_layout.id', name='plate_experiment_layout_fk', ondelete='cascade'))
    experiment_layout = relationship(
        ExperimentLayout,
        backref=backref('plates', uselist=True, cascade='delete,all'))
    barcode = Column(String(32), index=True)
    geid = Column(String(12), unique=True, nullable=False, index=True)
    description = Column(String)


guide_well_content_association = Table('guide_well_content_association', Base.metadata,
    Column('well_content_id', Integer, ForeignKey('well_content.id')),
    Column('guide_id', Integer, ForeignKey('guide.id'))
)


class WellContent(Base):
    __tablename__ = 'well_content'
    id = Column(Integer, primary_key=True) # is this id the replicate goup?
    clone_id = Column(Integer, ForeignKey('clone.id', name='well_content_clone_fk', ondelete='cascade'))
    clone = relationship(
        Clone,
        backref=backref('well_contents', uselist=True, cascade='delete,all'))
    guides = relationship('Guide', secondary=guide_well_content_association)
    replicate_group = Column(Integer, nullable=False, default=0)
    content_type = Column(Enum('WT', 'KO', 'BG', 'NM', 'SM', name='content_type'), nullable=False)
    is_control = Column(Boolean, nullable=False, default=False)


class Well(Base):
    __tablename__ = 'well'
    id = Column(Integer, primary_key=True)
    experiment_layout_id = Column(Integer, ForeignKey('experiment_layout.id', name='well_experiment_layout_fk', ondelete='cascade'), nullable=False)
    experiment_layout = relationship(
        ExperimentLayout,
        backref=backref('wells', uselist=True, cascade='delete,all'))
    well_content_id = Column(Integer, ForeignKey('well_content.id', name='well_well_content_fk', ondelete='cascade'), nullable=True)
    well_content = relationship(
        WellContent,
        backref=backref('wells', uselist=True, cascade='delete,all'))
    row = Column(String(1), nullable=False)
    column = Column(Integer, nullable=False)
    UniqueConstraint('experiment_layout_id', 'row', 'column', name='well_unique_in_layout')

    def isEmpty(self):
        return self.clone == None


class SequencingLibrary(Base):
    """
    template file: data/templates/YYYYMMDD_GEPXXXXX_SequencingLibrary.csv
    experiment_layout_geid	well_position	dna_source	slxid	library_type	barcode_size	sequencing_barcode	sequencing_sample_name
    """
    __tablename__ = 'sequencing_library'
    id = Column(Integer, primary_key=True)
    slxid = Column(String(8), nullable=False)
    library_type = Column(String(32))
    barcode_size = Column(Integer)


class SequencingLibraryContent(Base):
    __tablename__ = 'sequencing_library_content'
    id = Column(Integer, primary_key=True)
    well_id = Column(Integer, ForeignKey('well.id', name='well_sequencing_library_content_fk', ondelete='cascade'), nullable=True)
    well = relationship(
        Well,
        backref=backref('sequencing_library_contents', uselist=True, cascade='delete,all'))
    sequencing_library_id = Column(Integer, ForeignKey('sequencing_library.id', name='well_sequencing_library_fk', ondelete='cascade'), nullable=False)
    sequencing_library = relationship(
        SequencingLibrary,
        backref=backref('sequencing_library_contents', uselist=True, cascade='delete,all'))
    dna_source = Column(Enum('fixed cells', 'gDNA', 'non-fixed cells', name='dna_source'), nullable=False)
    sequencing_barcode = Column(String(20), nullable=False)
    sequencing_sample_name = Column(String(32), nullable=True)


class ProteinAbundance(Base):
    __tablename__ = 'abundance'
    id = Column(Integer, primary_key=True)
    well_id = Column(Integer, ForeignKey('well.id', name='abundance_well_fk', ondelete='cascade'), nullable=False)
    well = relationship(
        Well,
        backref=backref('abundances', uselist=True, cascade='delete,all'))
    plate_id = Column(Integer, ForeignKey('plate.id', name='abundance_plate_fk', ondelete='cascade'), nullable=False)
    plate = relationship(
        Plate,
        backref=backref('abundances', uselist=True, cascade='delete,all'))
    intensity_channel_700 = Column(Float)
    intensity_channel_800 = Column(Float)


class CellGrowth(Base):
    __tablename__ = 'growth'
    id = Column(Integer, primary_key=True)
    well_id = Column(Integer, ForeignKey('well.id', name='growth_well_fk', ondelete='cascade'), nullable=False)
    well = relationship(
        Well,
        backref=backref('growths', uselist=True, cascade='delete,all'))
    plate_id = Column(Integer, ForeignKey('plate.id', name='abundance_plate_fk', ondelete='cascade'), nullable=False)
    plate = relationship(
        Plate,
        backref=backref('growths', uselist=True, cascade='delete,all'))
    timestamp = Column(DateTime, nullable=False)
    hours = Column(Integer, nullable=False)
    confluence_percentage = Column(Float, nullable=False)
