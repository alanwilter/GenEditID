BASEDIR=$(dirname "$0")/..
cd $BASEDIR

source venv/bin/activate
echo 'Starting GenEditID WebApp'
pserve python/webapp/development.ini --reload
