"""
Python3 file of the genome-editing project
Created by Anne Pajon @pajanne on 08/03/2017
"""
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship
from sqlalchemy.orm import backref
from sqlalchemy import Table, Column, Integer, String, DateTime, Float, Boolean, Enum, ForeignKey, UniqueConstraint

Base = declarative_base()


class Project(Base):
    __tablename__ = 'project'
    id = Column(Integer, primary_key=True)
    name = Column(String(64), unique=True, nullable=False, index=True)
    description = Column(String(1024))


class ExperimentLayout(Base):
    __tablename__ = 'experiment_layout'
    id = Column(Integer, primary_key=True)
    project_id = Column(Integer, ForeignKey('project.id', name='experiment_layout_project_fk', ondelete='cascade'))
    project = relationship(
        Project,
        backref=backref('experiment_layouts', uselist=True, cascade='delete,all'))
    name = Column(String(32), index=True)
    description = Column(String)


class Target(Base):
    __tablename__ = 'target'
    id = Column(Integer, primary_key=True)
    project_id = Column(Integer, ForeignKey('project.id', name='target_project_fk', ondelete='cascade'))
    project = relationship(
        Project,
        backref=backref('targets', uselist=True, cascade='delete,all'))
    species = Column(String(32), nullable=False, index=True)
    chromosome = Column(String(32), nullable=False, index=True)
    name = Column(String(32), nullable=False, index=True)
    description = Column(String(1024))
    start = Column(Integer, nullable=False)
    end = Column(Integer, nullable=False)


class Guide(Base):
    __tablename__ = 'guide'
    id = Column(Integer, primary_key=True)
    target_id = Column(Integer, ForeignKey('target.id', name='guide_target_fk', ondelete='cascade'))
    target = relationship(
        Target,
        backref=backref('guides', uselist=True, cascade='delete,all'))
    amplicons = relationship("Amplicon", back_populates="guide")
    name = Column(String(32), nullable=False, index=True)
    description = Column(String(1024))
    location = Column(Integer, nullable=False)
    strand = Column(String(1), nullable=False)
    guide_sequence = Column(String(250), nullable=False)
    pam_sequence = Column(String(6), nullable=False)
    nuclease = Column(String(250), nullable=False)


class AmpliconSelection(Base):
    __tablename__ = 'amplicon_selection'
    guide_id = Column(Integer, ForeignKey('guide.id'), primary_key=True)
    guide = relationship("Guide", back_populates="amplicons")
    amplicon_id = Column(Integer, ForeignKey('amplicon.id'), primary_key=True)
    amplicon = relationship("Amplicon", back_populates="guides")
    name = Column(String(32), unique=True, nullable=False, index=True)
    description = Column(String(1024))
    score = Column(Integer)
    experiment_type = Column(String(32))
    donor_sequence = Column(String(250), nullable=True)


class Amplicon(Base):
    __tablename__ = 'amplicon'
    id = Column(Integer, primary_key=True)
    guides = relationship("Guide", back_populates="amplicon")
    name = Column(String(32), unique=True, nullable=False, index=True)
    description = Column(String(1024))
    target = Column(Enum('On', 'Off'))
    target_type = Column(Enum('gene', 'precursor', 'non-coding'))
    chromosome = Column(String(32), nullable=False, index=True)
    start = Column(Integer, nullable=False)
    end = Column(Integer, nullable=False)


class Primer(Base):
    __tablename__ = 'primer'
    id = Column(Integer, primary_key=True)
    amplicon_id = Column(Integer, ForeignKey('amplicon.id', name='primer_amplicon_fk', ondelete='cascade'))
    amplicon = relationship(
        Amplicon,
        backref=backref('primers', uselist=True, cascade='delete,all'))
    name = Column(String(32), nullable=False, index=True)
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
    control = Column(Boolean, nullable=False, default=False)


class Plate(Base):
    __tablename__ = 'plate'
    id = Column(Integer, primary_key=True)
    experiment_layout_id = Column(Integer, ForeignKey('experiment_layout.id', name='plate_experiment_layout_fk', ondelete='cascade'))
    experiment_layout = relationship(
        ExperimentLayout,
        backref=backref('plates', uselist=True, cascade='delete,all'))
    barcode = Column(String(32), index=True)
    description = Column(String)


guide_well_content_association = Table('guide_well_content_association', Base.metadata,
    Column('well_content_id', Integer, ForeignKey('well_content.id')),
    Column('guide_id', Integer, ForeignKey('guide.id'))
)


class WellContent(Base):
    __tablename__ = 'well_content'
    id = Column(Integer, primary_key=True) # is this id the replicate goup?
    replicate_group = Column(Integer, nullable=False, default=0)
    content_type = Column(Enum('wild_type', 'knockout', 'background', 'normaliser', 'empty', 'sample'), nullable=False)
    clone_id = Column(Integer, ForeignKey('clone.id', name='well_content_clone_fk', ondelete='cascade'))
    clone = relationship(
        Clone,
        backref=backref('well_contents', uselist=True, cascade='delete,all'))
    guides = relationship('Guide', secondary=guide_well_content_association)


class SequencingLibrary(Base):
    __tablename__ = 'sequencing_library'
    id = Column(Integer, primary_key=True)
    slxid = Column(String(8), nullable=False)
    library_type = Column(String(32))
    barcode_size = Column(Integer)


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
    sequencing_library_id = Column(Integer, ForeignKey('sequencing_library.id', name='well_sequencing_library_fk', ondelete='cascade'), nullable=True)
    sequencing_library = relationship(
        SequencingLibrary,
        backref=backref('wells', uselist=True, cascade='delete,all'))
    sequencing_barcode = Column(String(20))
    row = Column(String(1), nullable=False)
    column = Column(Integer, nullable=False)
    UniqueConstraint('experiment_layout_id', 'row', 'column', name='well_unique_in_layout')

    def isEmpty(self):
        return self.clone == None


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
