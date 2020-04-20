import os
import sys
import sqlalchemy

from geneditid.config import cfg
from geneditid.model import Base, Project

import log as logger


def main():

    log = logger.get_custom_logger(os.path.join(cfg['PROJECTS_FOLDER'], 'delete_project.log'))

    engine = sqlalchemy.create_engine(cfg['DATABASE_URI'])
    Base.metadata.bind = engine
    DBSession = sqlalchemy.orm.sessionmaker(bind=engine)
    session = DBSession()

    try:
        project = session.query(Project).filter(Project.geid == sys.argv[1]).first()
        session.delete(project)
        session.flush()
        session.commit()
        log.info('Project {} deleted'.format(sys.argv[1]))
    except Exception as e:
        log.exception(e)
        session.rollback()
    finally:
        session.close()


if __name__ == '__main__':
    main()
