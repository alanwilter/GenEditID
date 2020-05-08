from geneditid.config import yml_filepath, cfg

def test_yml_filepath():
    assert yml_filepath.endswith('geneditid.yml')

def test_cfg():
    assert 'DATABASE_URI' in cfg.keys()
    assert 'PROJECTS_FOLDER' in cfg.keys()
    assert 'FASTQ_SUBFOLDER' in cfg.keys()
