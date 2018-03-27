BASEDIR=$(dirname "$0")/..
cd $BASEDIR

source venv/bin/activate

for projectid in GEP00001 GEP00002 GEP00003 GEP00004 GEP00005 GEP00006 GEP00007 GEP00009 GEP00010
do
  echo 'Loading project' ${projectid}
  ./shell/load_project_${projectid}.sh
  echo '-----------------------------------------------------------------------'
done
