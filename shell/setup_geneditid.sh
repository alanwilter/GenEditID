BASEDIR=$(dirname "$0")/..
cd $BASEDIR

echo 'Creating Python3 virtual environment'
python3 -m venv venv
source venv/bin/activate

echo 'Installing GeneditID dependencies'
pip install -e python/.

echo 'Copying configuraton file'
cp python/dnascissors/geneditid.yml.sample python/dnascissors/geneditid.yml

echo 'Creating the database'
python python/scripts/create_db.py

echo 'Loading reference data'
python python/scripts/load_ref_data.py
