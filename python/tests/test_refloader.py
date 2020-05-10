from geneditid.config import cfg
from geneditid.connect import dbsession

from geneditid.loader import RefLoader
from geneditid.model import Genome
from geneditid.model import CellLine

class TestRefLoader:
    def test_load_genomes(self):
        ref_loader = RefLoader(dbsession)
        ref_loader.load_genomes()
        genomes = dbsession.query(Genome).all()
        for genome in genomes:
            found = False
            for genome_name in ref_loader.GENOMES:
                species = genome_name.split('[')[0][:-1]
                assembly = genome_name.split('[')[1][:-1]
                if genome.species == species and genome.assembly == assembly:
                    found = True
            assert found == True

    def test_load_celllines(self):
        ref_loader = RefLoader(dbsession)
        ref_loader.load_celllines()
        cell_lines = dbsession.query(CellLine).all()
        for cell_line in cell_lines:
            found = False
            for cell_line_name in ref_loader.CELL_LINES:
                if cell_line_name == cell_line.name:
                    found = True
            assert found == True
