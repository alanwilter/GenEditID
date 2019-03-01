BASEDIR=$(dirname "$0")/..
cd $BASEDIR

source venv/bin/activate
export PYTHONPATH=`pwd`/python

#shell/load_ref_data.sh

PROJNAME="GEP00013"
PROJDIR=data/${PROJNAME}

# echo 'Loading project' ${PROJNAME}
# python python/scripts/load_layout.py --layout=${PROJDIR}/${PROJNAME}.xlsx --geid=${PROJNAME}

for plateid in 01 02
do
  echo 'Loading data for plate' ${plateid}
  python python/scripts/load_cell_growth.py --plateid=${PROJNAME}_${plateid} --file=${PROJDIR}/${PROJNAME}_${plateid}_1in5_incu.txt
  #python python/scripts/load_cell_growth.py --plateid=${PROJNAME}_${plateid} --file=${PROJDIR}/${PROJNAME}_${plateid}_1in10_incu.txt
done

echo 'Loading sequencing variant results'
python python/scripts/load_variant_results.py --project_geid=${PROJNAME} --file=${PROJDIR}/${PROJNAME}.vardict.variants.xlsx --caller=VarDict
python python/scripts/load_variant_results.py --project_geid=${PROJNAME} --file=${PROJDIR}/${PROJNAME}.gatk.variants.xlsx --caller=HaplotypeCaller
python python/scripts/load_mutation_summary.py --project_geid=${PROJNAME}
