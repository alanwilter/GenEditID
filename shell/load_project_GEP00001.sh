BASEDIR=$(dirname "$0")/..
cd $BASEDIR

source venv/bin/activate
#export PYTHONPATH=`pwd`/python

shell/load_ref_data.sh

PROJNAME='GEP00001'
PROJDIR=data/${PROJNAME}

echo 'Loading project' ${PROJNAME}
python python/scripts/load_layout.py --layout=${PROJDIR}/${PROJNAME}.xlsx

for plateid in 01 02 03 04 05 06
do
  echo 'Loading data for plate' ${plateid}
  python python/scripts/load_protein_abundance.py --plateid=${PROJNAME}_${plateid}_ICW --file=${PROJDIR}/${PROJNAME}_${plateid}_ICW.csv
  python python/scripts/load_cell_growth.py --plateid=${PROJNAME}_${plateid}_incu --file=${PROJDIR}/${PROJNAME}_${plateid}_incu.txt
done
