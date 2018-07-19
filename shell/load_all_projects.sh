BASEDIR=$(dirname "$0")/..
cd $BASEDIR

source venv/bin/activate

for projectid in GEP00001 GEP00002 GEP00003 GEP00004 GEP00005 #GEP00006 GEP00007 GEP00009 GEP00010
do
  echo 'Loading project' ${projectid}
  ./shell/load_project_${projectid}.sh
  echo '-----------------------------------------------------------------------'
done

# Steps to complete the download of all current projects and associated data
# create project GEP00006 in webapp
# upload new project submission form and upload incucyte data on project GEP00006 edit page
# shell/load_project_GEP00006.sh # to load variant results
# shell/load_project_GEP00007.sh
# create project GEP00008 in webapp
# shell/load_project_GEP00009.sh
# shell/load_project_GEP00010.sh
