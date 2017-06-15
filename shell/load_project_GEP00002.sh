BASEDIR=$(dirname "$0")/..
cd $BASEDIR

source venv/bin/activate
#export PYTHONPATH=`pwd`/python

PROJDATE='20151101'
PROJNAME='GEP00002'
PROJDIR=data/${PROJDATE}_${PROJNAME}

echo 'Loading project' ${PROJNAME}
python python/scripts/load_layout.py --layout=${PROJDIR}/${PROJDATE}_${PROJNAME}.xlsx

for plateid in 01 02 03
do
  echo 'Loading data for plate' ${plateid}
  python python/scripts/load_protein_abundance.py --plateid=${PROJNAME}_${plateid}_ICW --file=${PROJDIR}/${PROJNAME}_${plateid}_ICW.csv
done

echo 'Loading sequencing variant results'
python python/scripts/load_variant_results.py --project_geid=${PROJNAME} --file=${PROJDIR}/SLX-13774.vardict.variants.xlsx --caller=VarDict
python python/scripts/load_variant_results.py --project_geid=${PROJNAME} --file=${PROJDIR}/SLX-13774.haplotypeCaller.variants.xlsx --caller=HaplotypeCaller
