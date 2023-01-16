import pytest
import sqlalchemy

from geneditid.config import cfg
from geneditid.connect import dbsession

from geneditid.loader import LoaderException
from geneditid.loader import ProjectLoader
from geneditid.loader import ProjectDataLoader

from geneditid.model import Project
from geneditid.model import Target
from geneditid.model import Guide
from geneditid.model import Amplicon
from geneditid.model import Primer

from geneditid.plotter import Plotter

class TestProjectDataLoader:
    def setup_method(self):
        project_loader = ProjectLoader(dbsession)
        project_loader.create_project('pytest')
        dbsession.commit()
        self.project = project_loader.project

    def teardown_method(self):
        dbsession.delete(self.project)
        dbsession.commit()

    def test_init_empty(self):
        loader = ProjectDataLoader(dbsession, self.project.geid, 'python/tests/pytest_empty.xlsx')
        assert loader.project.name == self.project.name
        assert loader.xls
        assert loader.xls.parse('Target').empty
        assert loader.xls.parse('Guide').empty
        assert loader.xls.parse('Amplicon').empty
        assert loader.xls.parse('Layout').empty

    def test_load_min(self):
        loader = ProjectDataLoader(dbsession, self.project.geid, 'python/tests/pytest_min.xlsx')
        assert not loader.xls.parse('Target').empty
        assert not loader.xls.parse('Guide').empty
        assert not loader.xls.parse('Amplicon').empty
        assert not loader.xls.parse('Layout').empty
        loader.load()
        dbsession.commit()
        target = dbsession.query(Target)\
                          .join(Project)\
                          .filter(Project.geid == self.project.geid)\
                          .filter(Target.name == 'test_target').one()
        assert target
        guide = dbsession.query(Guide)\
                         .join(Target)\
                         .join(Project)\
                         .filter(Project.geid == self.project.geid)\
                         .filter(Guide.name == 'test_guide')\
                         .one()
        assert guide
        amplicon = dbsession.query(Amplicon)\
                            .join(Project)\
                            .filter(Project.geid == self.project.geid)\
                            .filter(Amplicon.chromosome == 1)\
                            .filter(Amplicon.start == 1)\
                            .one()
        assert amplicon
        primer = dbsession.query(Primer)\
                          .join(Amplicon)\
                          .join(Project)\
                          .filter(Project.geid == self.project.geid)\
                          .filter(Amplicon.chromosome == 1)\
                          .filter(Amplicon.start == 1)\
                          .filter(Primer.sequence == 'AAAAA')\
                          .filter(Primer.strand == 'forward')\
                          .filter(Primer.start == 1)\
                          .filter(Primer.end == 10)\
                          .one()
        assert primer
        primer = dbsession.query(Primer)\
                          .join(Amplicon)\
                          .join(Project)\
                          .filter(Project.geid == self.project.geid)\
                          .filter(Amplicon.chromosome == 1)\
                          .filter(Amplicon.start == 1)\
                          .filter(Primer.sequence == 'CCCCC')\
                          .filter(Primer.strand == 'reverse')\
                          .filter(Primer.start == 101)\
                          .filter(Primer.end == 110)\
                          .one()
        assert primer

    def test_load_min_missingvalues(self):
        loader = ProjectDataLoader(dbsession, self.project.geid, 'python/tests/pytest_min_missingvalues.xlsx')
        with pytest.raises(LoaderException, match=r"Column target_gene_id needs a value in Target tab"):
            loader.load()

    def test_load_min_withlayouterror1(self):
        loader = ProjectDataLoader(dbsession, self.project.geid, 'python/tests/pytest_min_withlayouterror1.xlsx')
        with pytest.raises(sqlalchemy.exc.IntegrityError, match=r".*UNIQUE constraint failed: layout_content.row, layout_content.column, layout_content.layout_id.*"):
            try:
                loader.load()
            except sqlalchemy.exc.IntegrityError as e:
                dbsession.rollback()
                raise

    def test_load_min_withlayouterror2(self):
        loader = ProjectDataLoader(dbsession, self.project.geid, 'python/tests/pytest_min_withlayouterror2.xlsx')
        with pytest.raises(sqlalchemy.exc.IntegrityError, match=r".*UNIQUE constraint failed: layout_content.sequencing_barcode, layout_content.layout_id.*"):
            try:
                loader.load()
            except sqlalchemy.exc.IntegrityError as e:
                dbsession.rollback()
                raise

    def test_load_mouse(self):
        loader = ProjectDataLoader(dbsession, self.project.geid, 'python/tests/pytest_mouse.xlsx')
        assert not loader.xls.parse('Target').empty
        assert not loader.xls.parse('Guide').empty
        assert not loader.xls.parse('Amplicon').empty
        assert not loader.xls.parse('Layout').empty
        loader.load()
        dbsession.commit()
        target = dbsession.query(Target)\
                          .join(Project)\
                          .filter(Project.geid == self.project.geid)\
                          .filter(Target.name == 'test_mouse_target').one()
        assert target
        assert target.genome.name == 'Mus_musculus'

    @pytest.mark.parametrize(
        ("type_tuple", "input"),
        (
            # (('Mismatch', 'Substitution', 0.7),('5',66153585,'C','T')),
            # (('Insertion', 'FrameShift', 1.0),('5',66178781,'A','AA')))
            (('Mismatch', 'missense_variant', 0),('5',66153585,'C','T')),
            (('Insertion', 'frameshift_variant', 0),('5',66178781,'A','AA')),
            (('Deletion', 'missense_variant', 0),('16', 53704212,'CGAGAGCGCGAAGCTAA','C'))
            )
        )
    def test_get_variant_classification(self, type_tuple, input):
        p = Plotter(dbsession, self.project.geid)
        v = p.get_variant_classification(*input)
        assert type_tuple[1] in v
        assert v == type_tuple

    # def test_load_zebrafish(self):
    #     loader = ProjectDataLoader(dbsession, self.project.geid, 'python/tests/pytest_zebrafish.xlsx')
    #     assert not loader.xls.parse('Target').empty
    #     assert not loader.xls.parse('Guide').empty
    #     assert not loader.xls.parse('Amplicon').empty
    #     assert not loader.xls.parse('Layout').empty
    #     loader.load()
    #     dbsession.commit()
    #     target = dbsession.query(Target)\
    #                       .join(Project)\
    #                       .filter(Project.geid == self.project.geid)\
    #                       .filter(Target.name == 'test_mouse_target').one()
    #     assert target
    #     assert target.genome.name == 'Danio_rerio'
