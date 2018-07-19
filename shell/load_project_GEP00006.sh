BASEDIR=$(dirname "$0")/..
cd $BASEDIR

source venv/bin/activate
#export PYTHONPATH=`pwd`/python

shell/load_ref_data.sh

PROJNAME='GEP00006'
PROJDIR=data/${PROJNAME}

#echo 'Loading project' ${PROJNAME}
#python python/scripts/load_layout.py --layout=${PROJDIR}/${PROJNAME}.xlsx

#echo 'Loading data for plate 01'
#python python/scripts/load_cell_growth.py --plateid=${PROJNAME}_01_NGS --file=${PROJDIR}/JAG1_10.txt

echo 'Loading sequencing variant results'
python python/scripts/load_variant_results.py --project_geid=${PROJNAME} --file=${PROJDIR}/${PROJNAME}.vardict.variants.xlsx --caller=VarDict
python python/scripts/load_variant_results.py --project_geid=${PROJNAME} --file=${PROJDIR}/${PROJNAME}.gatk.variants.xlsx --caller=HaplotypeCaller
python python/scripts/load_mutation_summary.py --project_geid=${PROJNAME}
