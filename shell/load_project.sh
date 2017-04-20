BASEDIR=$(dirname "$0")/..
cd $BASEDIR

source venv/bin/activate
#export PYTHONPATH=`pwd`/python

PROJDATE='20170127'
PROJNAME='GEP00001'

echo 'Loading project' ${PROJNAME}
python python/scripts/load_layout.py --layout=data/${PROJDATE}_${PROJNAME}/${PROJDATE}_${PROJNAME}.xlsx

for plateid in 01 02 03 04 05 06
do
  echo 'Loading data for plate' ${plateid}
  python python/scripts/load_protein_abundance.py --plateid=${PROJNAME}_${plateid} --file=data/${PROJDATE}_${PROJNAME}/${PROJNAME}_${plateid}_ICW.csv
  python python/scripts/load_cell_growth.py --plateid=${PROJNAME}_${plateid} --file=data/${PROJDATE}_${PROJNAME}/${PROJNAME}_${plateid}_incu.txt
done
