"""
Python3 file of the genome-editing project
Created by Anne Pajon @pajanne on 08/03/2017
"""
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship
from sqlalchemy.orm import backref
from sqlalchemy import Column, Integer, String, DateTime, Float, Boolean, Enum, ForeignKey, UniqueConstraint

Base = declarative_base()


class Project(Base):
    __tablename__ = 'project'
    id = Column(Integer, primary_key=True)
    name = Column(String(64), unique=True, nullable=False, index=True)
    description = Column(String(1024))


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
    name = Column(String(32), nullable=False, index=True)
    description = Column(String(1024))
    location = Column(Integer, nullable=False)
    strand = Column(String(1), nullable=False)
    guide_sequence = Column(String(250), nullable=False)
    pam_sequence = Column(String(6), nullable=False)
    signature_protein = Column(String(250), nullable=False)


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


class Experiment(Base):
    __tablename__ = 'experiment'
    id = Column(Integer, primary_key=True)
    project_id = Column(Integer, ForeignKey('project.id', name='experiment_project_fk', ondelete='cascade'))
    project = relationship(
        Project,
        backref=backref('experiments', uselist=True, cascade='delete,all'))
    name = Column(String(32), index=True)
    description = Column(String)


class Plate(Base):
    __tablename__ = 'plate'
    id = Column(Integer, primary_key=True)
    experiment_id = Column(Integer, ForeignKey('experiment.id', name='plate_experiment_fk', ondelete='cascade'))
    experiment = relationship(
        Experiment,
        backref=backref('plates', uselist=True, cascade='delete,all'))
    name = Column(String(32), index=True)
    description = Column(String)


class Well(Base):
    __tablename__ = 'well'
    id = Column(Integer, primary_key=True)
    plate_id = Column(Integer, ForeignKey('plate.id', name='well_plate_fk', ondelete='cascade'), nullable=False)
    plate = relationship(
        Plate,
        backref=backref('wells', uselist=True, cascade='delete,all'))
    clone_id = Column(Integer, ForeignKey('clone.id', name='well_clone_fk', ondelete='cascade'))
    clone = relationship(
        Clone,
        backref=backref('wells', uselist=True, cascade='delete,all'))
    target_id = Column(Integer, ForeignKey('target.id', name='well_target_fk', ondelete='cascade'))
    target = relationship(
        Target,
        backref=backref('wells', uselist=True, cascade='delete,all'))
    guide_id = Column(Integer, ForeignKey('guide.id', name='well_guide_fk', ondelete='cascade'))
    guide = relationship(
        Guide,
        backref=backref('wells', uselist=True, cascade='delete,all'))
    row = Column(String(1), nullable=False)
    column = Column(Integer, nullable=False)
    replicate_group = Column(Integer, nullable=False, default=0)
    well_type = Column(Enum('wild_type', 'knockout', 'background', 'normaliser', 'empty', 'sample'), nullable=False)
    UniqueConstraint('plate_id', 'row', 'column', name='well_unique_in_plate')

    def isEmpty(self):
        return self.clone == None


class ProteinAbundance(Base):
    __tablename__ = 'abundance'
    id = Column(Integer, primary_key=True)
    well_id = Column(Integer, ForeignKey('well.id', name='abundance_well_fk', ondelete='cascade'), nullable=False)
    well = relationship(
        Well,
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
    timestamp = Column(DateTime, nullable=False)
    hours = Column(Integer, nullable=False)
    confluence_percentage = Column(Float, nullable=False)
