BASEDIR=$(dirname "$0")/..
cd $BASEDIR

source venv/bin/activate
#export PYTHONPATH=`pwd`/python

echo 'Loading reference data'
python python/scripts/load_ref_data.py
