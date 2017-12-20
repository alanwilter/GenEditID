BASEDIR=$(dirname "$0")/..
cd $BASEDIR

source venv/bin/activate
#export PYTHONPATH=`pwd`/python

shell/load_ref_data.sh

PROJNAME='GEP00005'
PROJDIR=data/${PROJNAME}

echo 'Loading project' ${PROJNAME}
python python/scripts/load_layout.py --layout=${PROJDIR}/${PROJNAME}.xlsx

#echo 'Loading sequencing variant results'
#python python/scripts/load_variant_results.py --project_geid=${PROJNAME} --file=${PROJDIR}/${PROJNAME}.vardict.variants.xlsx --caller=VarDict
#python python/scripts/load_variant_results.py --project_geid=${PROJNAME} --file=${PROJDIR}/${PROJNAME}.gatk.variants.xlsx --caller=HaplotypeCaller
#python python/scripts/load_mutation_summary.py --project_geid=${PROJNAME}

# Loading sequencing variant results
# dnascissors              ERROR   : There is no sequencing library content for MixA_KO_A1 sample with FLD0001 barcode
# Traceback (most recent call last):
#   File "python/scripts/load_variant_results.py", line 32, in main
#     loader.load()
#   File "/Users/pajon01/workspace/genome-editing/python/dnascissors/loader.py", line 729, in load
#     self.load_sheet('SNVs', 'SNV')
#   File "/Users/pajon01/workspace/genome-editing/python/dnascissors/loader.py", line 695, in load_sheet
#     raise LoaderException("There is no sequencing library content for {} sample with {} barcode".format(row.sample, row.barcode))
