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
    name = Column(String, unique = True, nullable = False)


class Target(Base):
    __tablename__ = 'target'
    id = Column(Integer, primary_key=True)
    project_id = Column(Integer, ForeignKey('project.id', name = 'target_project_fk', ondelete = 'cascade'))
    # Propagate the deletion of a Project onto its targets
    project = relationship(
        Project,
        backref=backref('targets', uselist=True, cascade='delete,all'))
    name = Column(String, nullable = False)
    description = Column(String)


class Oligo(Base):
    __tablename__ = 'oligo'
    id = Column(Integer, primary_key=True)
    target_id = Column(Integer, ForeignKey('target.id', name = 'oligo_target_fk', ondelete = 'cascade'))
    target = relationship(
        Target,
        backref=backref('oligos', uselist=True, cascade='delete,all'))
    name = Column(String, nullable = False)
    description = Column(String)


class CellLine(Base):
    __tablename__ = 'cell_line'
    id = Column(Integer, primary_key=True)
    name = Column(String, unique=True)


class Clone(Base):
    __tablename__ = 'clone'
    id = Column(Integer, primary_key=True)
    oligo_id = Column(Integer, ForeignKey('oligo.id', name='clone_oligo_fk', ondelete = 'cascade'))
    oligo = relationship(
        Oligo,
        backref=backref('clones', uselist=True, cascade='delete,all'))
    cell_line_id = Column(Integer, ForeignKey('cell_line.id', name='clone_cellline_fk', ondelete = 'cascade'))
    cell_line = relationship(
        CellLine,
        backref=backref('clones', uselist=True, cascade='delete,all'))
    name = Column(String, nullable = False, unique = True)
    description = Column(String)

class Plate(Base):
    __tablename__ = 'plate'
    id = Column(Integer, primary_key=True)
    name = Column(String)
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
    location = Column(String)
    empty = Column(Boolean)
    UniqueConstraint('plate_id', 'location', name='well_unique_in_plate')

"""
class Measurement(Base):
    __tablename__ = 'measurement'
    id = Column(Integer, primary_key=True)
    well_id = Column(Integer, ForeignKey('well.id', name = 'measurement_well_fk', ondelete = 'cascade'))
    well = relationship(
        Well,
        backref=backref('measurements', uselist=True, cascade='delete,all'))
    measurement_type = Column(String)
    unit = Column(String)
    date = Column(DateTime)
    value = Column(Float)
"""
