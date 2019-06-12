BASEDIR=$(dirname "$0")/..
cd $BASEDIR

source venv/bin/activate
echo 'Start GeneditID WebApp'
pserve python/webapp/development.ini --reload
