BASEDIR=$(dirname "$0")/..
cd $BASEDIR

source venv/bin/activate
#export PYTHONPATH=`pwd`/python

PROJDATE='20170201'
PROJNAME='GEP00003'
PROJDIR=data/${PROJDATE}_${PROJNAME}

echo 'Loading project' ${PROJNAME}
python python/scripts/load_layout.py --layout=${PROJDIR}/${PROJDATE}_${PROJNAME}.xlsx
