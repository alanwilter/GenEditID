from geneditid.connect import engine, dbsession

def test_engine():
    assert engine

def test_dbsession():
    assert dbsession
