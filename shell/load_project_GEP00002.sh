BASEDIR=$(dirname "$0")/..
cd $BASEDIR

source venv/bin/activate
#export PYTHONPATH=`pwd`/python

shell/load_ref_data.sh

PROJNAME='GEP00002'
PROJDIR=data/${PROJNAME}

echo 'Loading project' ${PROJNAME}
python python/scripts/load_layout.py --layout=${PROJDIR}/${PROJNAME}.xlsx

for plateid in 01 02 03
do
  echo 'Loading data for plate' ${plateid}
  python python/scripts/load_protein_abundance.py --plateid=${PROJNAME}_${plateid}_ICW --file=${PROJDIR}/${PROJNAME}_${plateid}_ICW.csv
done

echo 'Loading sequencing variant results'
python python/scripts/load_variant_results.py --project_geid=${PROJNAME} --file=${PROJDIR}/SLX-13774.vardict.variants.xlsx --caller=VarDict
python python/scripts/load_variant_results.py --project_geid=${PROJNAME} --file=${PROJDIR}/SLX-13774.haplotypeCaller.variants.xlsx --caller=HaplotypeCaller
python python/scripts/load_mutation_summary.py --project_geid=${PROJNAME}
