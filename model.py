"""
Python3 file of the genome-editing project
Created by Anne Pajon @pajanne on 08/03/2017
"""
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship
from sqlalchemy.orm import backref
from sqlalchemy import Column
from sqlalchemy import Integer
from sqlalchemy import String
from sqlalchemy import DateTime

Base = declarative_base()


class Project(Base):
    __tablename__ = 'project'
    id = Column(Integer, primary_key=True)
    name = Column(String, unique=True)


class Target(Base):
    __tablename__ = 'target'
    id = Column(Integer, primary_key=True)
    # Propagate the deletion of a Project onto its targets
    project = relationship(
        Project,
        backref=backref('targets', uselist=True, cascade='delete,all'))
    name = Column(String)
    description = Column(String)


class Oligo(Base):
    __tablename__ = 'oligo'
    id = Column(Integer, primary_key=True)
    target = relationship(
        Target,
        backref=backref('oligos', uselist=True, cascade='delete,all'))
    name = Column(String)
    description = Column(String)


class CellLine(Base):
    __tablename__ = 'cell_line'
    id = Column(Integer, primary_key=True)
    name = Column(String, unique=True)


class Clone(Base):
    __tablename__ = 'clone'
    id = Column(Integer, primary_key=True)
    oligo = relationship(
        Oligo,
        backref=backref('clones', uselist=True, cascade='delete,all'))
    cell_line = relationship(
        CellLine,
        backref=backref('clones', uselist=True, cascade='delete,all'))
    name = Column(String)
    description = Column(String)


class Plate(Base):
    __tablename__ = 'plate'
    id = Column(Integer, primary_key=True)
    name = Column(String)
    description = Column(String)


class PlateLayout(Base):
    __tablename__ = 'plate_layout'
    id = Column(Integer, primary_key=True)
    plate = relationship(
        Plate,
        backref=backref('plate_layouts', uselist=True, cascade='delete,all'))
    clone = relationship(
        Clone,
        backref=backref('plate_layouts', uselist=True, cascade='delete,all'))
    location = Column(String)


class Measurement(Base):
    __tablename__ = 'measurement'
    id = Column(Integer, primary_key=True)
    plate_layout = relationship(
        PlateLayout,
        backref=backref('measurements', uselist=True, cascade='delete,all'))
    measurement_type = Column(String)
    unit = Column(String)
    date = Column(DateTime)
