#BASEDIR=$(dirname "$0")/..
BASEDIR="$( cd "$(dirname "$0")" ; pwd -P )"/..
cd $BASEDIR

echo 'Creating Python3 virtual environment'
python3 -m venv venv
source venv/bin/activate

echo 'Installing GenEditID dependencies'
pip install -e python/.

echo 'Copying GenEditID configuraton file'
cp python/geneditid/geneditid.yml.sample python/geneditid/geneditid.yml

echo 'Creating the GenEditID database and the projects folder'
python python/geneditidtools/create_db.py

echo 'Loading GenEditID reference data'
python python/geneditidtools/load_ref_data.py

echo 'Installing Homo Sapiens reference genome'
pyensembl install --release 95 --species homo_sapiens

cd python/geneditidapp/static/
if [ ! -h GEPXXXXX.xlsx ]; then
  echo 'Create symlink to template file for accessing it from the WebApp'
  ln -s ../../../data/templates/GEPXXXXX.xlsx
fi
cd $BASEDIR
