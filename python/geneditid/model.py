import os
import glob
from geneditid.config import cfg

from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship
from sqlalchemy.orm import backref
from sqlalchemy import Table, Column, Integer, String, Date, DateTime, Float, Boolean, Enum, ForeignKey, UniqueConstraint

Base = declarative_base()


class Genome(Base):
    # Ref from https://www.ncbi.nlm.nih.gov/grc
    __tablename__ = 'genome'
    __table_args__ = ( UniqueConstraint('species', 'assembly', name='unique_genome_species_assembly'),
                     )
    id = Column(Integer, primary_key=True)
    species = Column(String(32), nullable=False, index=True)
    assembly = Column(String(32), nullable=False, index=True)

    @property
    def fa_file(self):
        return os.path.join(cfg['DATA_FOLDER'], cfg['REF_SUBFOLDER'], '{}.{}.dna.toplevel.fa.gz'.format(self.species.replace(' ', '_'), self.assembly))


class CellLine(Base):
    # Ref from http://web.expasy.org/cellosaurus/
    # Data ftp://ftp.expasy.org/databases/cellosaurus/cellosaurus.txt
    __tablename__ = 'cell_line'
    __table_args__ = ( UniqueConstraint('name', name='cell_line_name'),
                     )
    id = Column(Integer, primary_key=True)
    name = Column(String(32), nullable=False, index=True)


class Project(Base):
    __tablename__ = 'project'
    __table_args__ = ( UniqueConstraint('geid', name='unique_project_geid'),
                       UniqueConstraint('name', name='unique_project_name')
                     )
    id = Column(Integer, primary_key=True)
    geid = Column(String(8), nullable=False, index=True)
    name = Column(String(64), nullable=False, index=True)
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
    def project_folder(self):
        return os.path.join(cfg['PROJECTS_FOLDER'], self.geid)

    @property
    def is_abundance_data_available(self):
        if self.layouts:
            for layout in self.layouts:
                if layout.layout_contents:
                    for layout_content in layout.layout_contents:
                        if len(layout_content.abundances) > 0:
                            return True
        return False

    @property
    def is_growth_data_available(self):
        if self.layouts:
            for layout in self.layouts:
                if layout.layout_contents:
                    for layout_content in layout.layout_contents:
                        if len(layout_content.growths) > 0:
                            return True
        return False

    @property
    def is_sequencing_data_available(self):
        if glob.glob(os.path.join(self.project_folder, cfg['FASTQ_SUBFOLDER'], '*.fq.gz')):
            # TODO check that all samples have a fastq files on disk
            return True
        return False


class Target(Base):
    __tablename__ = 'target'
    __table_args__ = ( UniqueConstraint('name', 'genome_id', 'project_id', name='unique_target_in_project'),
                     )
    id = Column(Integer, primary_key=True)
    project_id = Column(Integer, ForeignKey('project.id', name='target_project_fk', ondelete='CASCADE'), nullable=False)
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


class Guide(Base):
    __tablename__ = 'guide'
    __table_args__ = ( UniqueConstraint('name', 'genome_id', 'target_id', 'project_id', name='unique_guide_in_target_within_project'),
                     )
    id = Column(Integer, primary_key=True)
    project_id = Column(Integer, ForeignKey('project.id', name='guide_project_fk', ondelete='CASCADE'), nullable=False)
    project = relationship(
        Project,
        backref=backref('guides', uselist=True, cascade='delete,all'))
    target_id = Column(Integer, ForeignKey('target.id', name='guide_target_fk', ondelete='CASCADE'), nullable=False)
    target = relationship(
        Target,
        backref=backref('guides', uselist=True, cascade='delete,all'))
    genome_id = Column(Integer, ForeignKey('genome.id', name='guide_genome_fk'), nullable=False)
    genome = relationship(Genome)
    name = Column(String(32), nullable=False, index=True)
    guide_sequence = Column(String(250), nullable=False)
    pam_sequence = Column(String(6), nullable=True)
    activity = Column(Integer, nullable=True)
    exon = Column(Integer, nullable=True)
    nuclease = Column(String(250), nullable=True)


class GuideMismatch(Base):
    __tablename__ = 'guide_mismatch'
    id = Column(Integer, primary_key=True)
    guide_id = Column(Integer, ForeignKey('guide.id', name='guide_mismatch_guide_fk', ondelete='CASCADE'), nullable=False)
    guide = relationship(
        Guide,
        backref=backref('guide_mismatches', uselist=True, cascade='delete,all'))
    is_off_target_coding_region = Column(Boolean, nullable=False)
    number_of_mismatches = Column(Integer, nullable=False)
    number_of_off_targets = Column(Integer, nullable=False)


class Amplicon(Base):
    __tablename__ = 'amplicon'
    id = Column(Integer, primary_key=True)
    project_id = Column(Integer, ForeignKey('project.id', name='amplicon_project_fk', ondelete='CASCADE'), nullable=False)
    project = relationship(
        Project,
        backref=backref('amplicons', uselist=True, cascade='delete,all'))
    genome_id = Column(Integer, ForeignKey('genome.id', name='amplicon_genome_fk'), nullable=False)
    genome = relationship(Genome)
    guide_id = Column(Integer, ForeignKey('guide.id', name="amplicon_guide_fk", ondelete='CASCADE'), nullable=False)
    guide = relationship(
        Guide,
        backref=backref('amplicons', uselist=True, cascade='delete,all'))
    dna_feature = Column(Enum('gene', 'precursor', 'non-coding', name='dna_feature'))
    chromosome = Column(String(32), nullable=False, index=True)
    start = Column(Integer, nullable=False)
    end = Column(Integer, nullable=False)
    experiment_type = Column(String(32), nullable=False)
    guide_location = Column(Integer, nullable=False)
    guide_strand = Column(Enum('forward', 'reverse', name='strand'), nullable=False)
    is_on_target = Column(Boolean, nullable=False)
    score = Column(Integer)
    description = Column(String(1024))

    @property
    def name(self):
        return '{}_chr{}_{}'.format(self.guide.target.name, self.chromosome, self.start)

    @property
    def coordinates(self):
        return 'chr{}:{}-{}'.format(self.chromosome, self.start, self.end)

    @property
    def strand(self):
        if self.guide_strand == 'reverse':
            return '-'
        return '+'

    @property
    def fprimer(self):
        for primer in self.primers:
            if primer.strand == 'forward':
                return primer

    @property
    def rprimer(self):
        for primer in self.primers:
            if primer.strand == 'reverse':
                return primer

class Primer(Base):
    __tablename__ = 'primer'
    id = Column(Integer, primary_key=True)
    amplicon_id = Column(Integer, ForeignKey('amplicon.id', name='primer_amplicon_fk', ondelete='CASCADE'), nullable=False)
    amplicon = relationship(
        Amplicon,
        backref=backref('primers', uselist=True, cascade='delete,all'))
    genome_id = Column(Integer, ForeignKey('genome.id', name='primer_genome_fk'), nullable=False)
    genome = relationship(Genome)
    sequence = Column(String(250), nullable=False)
    strand = Column(Enum('forward', 'reverse', name='strand'), nullable=False)
    start = Column(Integer, nullable=False)
    end = Column(Integer, nullable=False)


class Clone(Base):
    __tablename__ = 'clone'
    __table_args__ = ( UniqueConstraint('name', 'project_id', name='unique_clone_in_project'),
                     )
    id = Column(Integer, primary_key=True)
    project_id = Column(Integer, ForeignKey('project.id', name='clone_project_fk', ondelete='CASCADE'), nullable=False)
    project = relationship(
        Project,
        backref=backref('clones', uselist=True, cascade='delete,all'))
    cell_line_id = Column(Integer, ForeignKey('cell_line.id', name='clone_cell_line_fk'))
    cell_line = relationship(CellLine)
    name = Column(String(32), nullable=False, index=True)
    cell_pool = Column(String(32))
    description = Column(String(1024))


class Layout(Base):
    __tablename__ = 'layout'
    __table_args__ = ( UniqueConstraint('geid', 'project_id', name='unique_layout_in_project'),
                     )
    id = Column(Integer, primary_key=True)
    project_id = Column(Integer, ForeignKey('project.id', name='layout_project_fk', ondelete='CASCADE'), nullable=False)
    project = relationship(
        Project,
        backref=backref('layouts', uselist=True, cascade='delete,all'))
    geid = Column(String(12), nullable=False, index=True)


class LayoutContent(Base):
    __tablename__ = 'layout_content'
    __table_args__ = ( UniqueConstraint('row', 'column', 'layout_id', name='layout_content_location_unique_in_layout'),
                       UniqueConstraint('sequencing_barcode', 'layout_id', name='sequencing_barcode_unique_in_layout'),
                     )
    id = Column(Integer, primary_key=True)
    layout_id = Column(Integer, ForeignKey('layout.id', name='layout_content_layout_fk', ondelete='CASCADE'), nullable=False)
    layout = relationship(
        Layout,
        backref=backref('layout_contents', uselist=True, cascade='delete,all'))
    guide_id = Column(Integer, ForeignKey('guide.id', name="layout_content_guide_fk", ondelete='CASCADE'))
    guide = relationship(
        Guide,
        backref=backref('layout_contents', uselist=True, cascade='delete,all'))
    clone_id = Column(Integer, ForeignKey('clone.id', name='layout_content_clone_fk', ondelete='CASCADE'))
    clone = relationship(
        Clone,
        backref=backref('layout_contents', uselist=True, cascade='delete,all'))
    row = Column(String(1), nullable=False)
    column = Column(Integer, nullable=False)
    sequencing_project_id = Column(String(12))
    sequencing_library_type = Column(String(32))
    sequencing_barcode = Column(String(20))
    sequencing_dna_source = Column(Enum('fixed cells', 'gDNA', 'non-fixed cells', 'water', name='dna_source'))
    sequencing_sample_name = Column(String(32))
    content_type = Column(Enum('wild-type', 'knock-out', 'background', 'normalisation', 'sample', 'empty-vector', 'empty', name='content_type'), nullable=False, default='sample')
    is_control = Column(Boolean, nullable=False, default=False)
    replicate_group = Column(Integer, nullable=False, default=0)

    def is_empty(self):
        return (self.content_type == 'empty' and self.guide == None)

    @property
    def position(self):
        return '{}{}'.format(self.row, self.column)


class Plate(Base):
    __tablename__ = 'plate'
    __table_args__ = ( UniqueConstraint('barcode', name='unique_plate_barcode'),
                       UniqueConstraint('name', 'layout_id', name='unique_plate_name_in_layout')
                     )
    id = Column(Integer, primary_key=True)
    layout_id = Column(Integer, ForeignKey('layout.id', name='plate_layout_fk', ondelete='CASCADE'), nullable=False)
    layout = relationship(
        Layout,
        backref=backref('plates', uselist=True, cascade='delete,all'))
    name = Column(String(32), nullable=False, index=True)
    barcode = Column(String(32), index=True)
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
    def plate_type(self):
        types = {'abundance': self.is_abundance_plate,
                 'growth': self.is_growth_plate,
                 }
        return ','.join([k for k, v in types.items() if v is True])


class ProteinAbundance(Base):
    __tablename__ = 'abundance'
    id = Column(Integer, primary_key=True)
    layout_content_id = Column(Integer, ForeignKey('layout_content.id', name='abundance_layout_content_fk', ondelete='CASCADE'), nullable=False)
    layout_content = relationship(
        LayoutContent,
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
    layout_content_id = Column(Integer, ForeignKey('layout_content.id', name='growth_layout_content_fk', ondelete='CASCADE'), nullable=False)
    layout_content = relationship(
        LayoutContent,
        backref=backref('growths', uselist=True, cascade='delete,all'))
    plate_id = Column(Integer, ForeignKey('plate.id', name='abundance_plate_fk', ondelete='CASCADE'), nullable=False)
    plate = relationship(
        Plate,
        backref=backref('growths', uselist=True, cascade='delete,all'))
    timestamp = Column(DateTime, nullable=False)
    hours = Column(Integer, nullable=False)
    confluence_percentage = Column(Float, nullable=False)
