"""
Python3 file of the genome-editing project
Created by Anne Pajon @pajanne on 08/03/2017
"""
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship
from sqlalchemy.orm import backref
from sqlalchemy import Column, Integer, String, DateTime, Float, Boolean, ForeignKey, UniqueConstraint

Base = declarative_base()


class Project(Base):
    __tablename__ = 'project'
    id = Column(Integer, primary_key=True)
    name = Column(String(64), unique = True, nullable = False, index = True)


class Target(Base):
    __tablename__ = 'target'
    id = Column(Integer, primary_key=True)
    project_id = Column(Integer, ForeignKey('project.id', name = 'target_project_fk', ondelete = 'cascade'))
    # Propagate the deletion of a Project onto its targets
    project = relationship(
        Project,
        backref=backref('targets', uselist=True, cascade='delete,all'))
    name = Column(String(32), nullable = False, index = True)
    description = Column(String(1024))
    start = Column(Integer, nullable = False)
    end = Column(Integer, nullable = False)


class SGRNA(Base):
    __tablename__ = 'sgrna'
    id = Column(Integer, primary_key=True)
    target_id = Column(Integer, ForeignKey('target.id', name = 'sgrna_target_fk', ondelete = 'cascade'))
    target = relationship(
        Target,
        backref=backref('sgrnas', uselist=True, cascade='delete,all'))
    name = Column(String(32), nullable = False, index = True)
    description = Column(String(1024))
    location = Column(Integer, nullable = False)


class CellLine(Base):
    __tablename__ = 'cell_line'
    id = Column(Integer, primary_key = True)
    name = Column(String(32), unique = True, nullable = False, index = True)


class Clone(Base):
    __tablename__ = 'clone'
    id = Column(Integer, primary_key=True)
    sgrna_id = Column(Integer, ForeignKey('sgrna.id', name='clone_sgrna_fk', ondelete = 'cascade'))
    sgrna = relationship(
        SGRNA,
        backref=backref('clones', uselist=True, cascade='delete,all'))
    cell_line_id = Column(Integer, ForeignKey('cell_line.id', name='clone_cellline_fk', ondelete = 'cascade'))
    cell_line = relationship(
        CellLine,
        backref=backref('clones', uselist=True, cascade='delete,all'))
    name = Column(String(32), nullable = False, unique = True, index = True)
    description = Column(String(1024))
    control = Column(Boolean, nullable = False, default = False)


class Plate(Base):
    __tablename__ = 'plate'
    id = Column(Integer, primary_key=True)
    name = Column(String(32), index = True)
    description = Column(String)


class Well(Base):
    __tablename__ = 'well'
    id = Column(Integer, primary_key=True)
    plate_id = Column(Integer, ForeignKey('plate.id', name='well_plate_fk', ondelete = 'cascade'), nullable = False)
    plate = relationship(
        Plate,
        backref=backref('wells', uselist=True, cascade='delete,all'))
    clone_id = Column(Integer, ForeignKey('clone.id', name='well_clone_fk', ondelete = 'cascade'))
    clone = relationship(
        Clone,
        backref=backref('wells', uselist=True, cascade='delete,all'))
    row = Column(String(1), nullable = False)
    column = Column(Integer, nullable = False)
    icw_700 = Column(Float)
    icw_800 = Column(Float)
    UniqueConstraint('plate_id', 'row', 'column', name='well_unique_in_plate')
    
    def isEmpty(self):
        return self.clone == None


class IncucyteGrowth(Base):
    __tablename__ = 'incucyte'
    id = Column(Integer, primary_key=True)
    well_id = Column(Integer, ForeignKey('well.id', name='incucyte_well_fk', ondelete = 'cascade'), nullable = False)
    well = relationship(
        Well,
        backref=backref('incucyte', uselist=True, cascade='delete,all'))
    timestamp = Column(DateTime, nullable = False)
    hours = Column(Integer, nullable = False)
    phased_object_confluence = Column(Float, nullable = False)
