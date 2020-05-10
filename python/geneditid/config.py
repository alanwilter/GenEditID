import os
import yaml

env = os.environ.get('PYENV')

if env == 'test':
    yml_filepath = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'geneditid_test.yml')
else:
    yml_filepath = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'geneditid.yml')

with open(yml_filepath, 'r') as yml_file:
    cfg = yaml.load(yml_file, Loader=yaml.FullLoader)
    if 'sqlite' in cfg['DATABASE_URI']:
        db_path = os.path.join(os.path.dirname(__file__), '..', '..', cfg['DATABASE_URI'].split('/')[3])
        cfg['DATABASE_URI'] = 'sqlite:///{}'.format(db_path)

# create projects folder
if not os.path.exists(cfg['PROJECTS_FOLDER']):
    os.makedirs(cfg['PROJECTS_FOLDER'])
