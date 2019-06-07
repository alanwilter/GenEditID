# GenEditID WebApp

## Getting started

- Change directory into your newly created project.
  ```
  git clone https://github.com/GenEditID/editID.git
  cd editID/python/webapp
  ```
- Create a Python virtual environment.
  ```
  python3 -m venv wenv
  ```
- Activate the virtual environment for the webapp
  ```
  source wenv/bin/activate
  ```
- Upgrade packaging tools.
  ```
  pip install --upgrade pip
  ```
- Install dependencies
  ```
  pip install -r requirements-webapp.txt
  ```
- Setting up a local postgres database
  ```
  createdb geneditid
  createuser gene
  psql geneditid
  geneediting=# GRANT CREATE ON SCHEMA public TO "gene";
  geneediting=# GRANT USAGE ON SCHEMA public TO "gene";
  geneediting=# ALTER USER gene WITH PASSWORD 'gene';
  # edit `crispr.yml` to set DATABASE_URI: "postgresql://gene:gene@localhost/geneditid"
  python python/scripts/create_db.py
  ```
- Run your webapp
  ```
  pserve development.ini --reload
  ```
  Go to http://localhost:8080
- Create a new project and load the "data and experiment layout submission" excel (.xls) file


## Running the WebApp after first setup

Edit `python/webapp/development.ini` to check you are pointing to the right database

```
cd editID/python/webapp
source wenv/bin/activate
pserve development.ini --reload
```

Go to http://localhost:8080


## Setting up a linux server for production

https://docs.pylonsproject.org/projects/pyramid/en/latest/tutorials/modwsgi/index.html

Once you have Apache installed, install mod_wsgi.
https://github.com/GrahamDumpleton/mod_wsgi

```bash
sudo yum install httpd httpd-devel
sudo yum install python34-devel

sudo su - ge
source venv/bin/activate
pip install mod_wsgi
mod_wsgi-express start-server
```

Verify that the installation worked by pointing a browser to the expected url.

Create `editID/python/webapp/pyramid.wsgi` with this content
```python
from pyramid.paster import get_app, setup_logging
ini_path = '/home/ge/editID/python/webapp/production.ini'
setup_logging(ini_path)
application = get_app(ini_path, 'main')
```

Configure Apache, with file `vim /etc/httpd/conf.d/geneditid.conf`
NB. remove specific configuration for the time being otherwise links generated
by Pyramid become absolute and come up as http://localhost:8080/ which is a problem.
Restart Apache
```bash
sudo /sbin/service httpd restart
```

Starting the web app:
```
sudo su - ge
source venv/bin/activate
cd editID/python/webapp
mod_wsgi-express start-server pyramid.wsgi --port 8080 --user ge --group users --server-root=/home/ge/
```

Restarting the web app:

```bash
sudo su - ge
source venv/bin/activate

# kill the current processes
ps -eaf | grep wsgi

# update app
cd genome-editing
git pull

# update dependencies
cd python/webapp
pip install -r requirements-webapp.txt

# start
./prod_startwebapp &
```

Log file `/home/ge/error_log`

NB. currently running on port 8080, apache config needs to be changed.
