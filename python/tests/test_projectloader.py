from geneditid.config import cfg
from geneditid.connect import dbsession

from geneditid.model import Project
from geneditid.loader import ProjectLoader

class TestProjectLoader:
    def test_set_next_project_geid(self):
        project_loader = ProjectLoader(dbsession)
        project_loader.set_next_project_geid()
        last_project = dbsession.query(Project).order_by(Project.id.desc()).first()
        last_project_geid = 'GEP00000'
        if last_project:
            last_project_geid = last_project.geid
        assert int(last_project_geid[3:]) + 1 == int(project_loader.project_geid[3:])

    def test_create_project(self):
        project_loader = ProjectLoader(dbsession)
        project_loader.create_project('pytest project', 'knock-out')
        dbsession.commit()
        last_project = dbsession.query(Project).order_by(Project.id.desc()).first()
        assert last_project.name == project_loader.project.name
        project_loader.delete_project(project_loader.project.geid)
        dbsession.commit()
        deleted_project = dbsession.query(Project).filter(Project.geid == project_loader.project.geid).first()
        assert not deleted_project

    def test_delete_project(self):
        project_loader = ProjectLoader(dbsession)
        project_loader.create_project('pytest project', 'knock-out')
        project_geid_todelete = project_loader.project_geid
        dbsession.commit()
        last_project = dbsession.query(Project).order_by(Project.id.desc()).first()
        assert project_geid_todelete == last_project.geid
        project_loader.delete_project(project_geid_todelete)
        dbsession.commit()
        deleted_project = dbsession.query(Project).filter(Project.geid == project_geid_todelete).first()
        assert not deleted_project
