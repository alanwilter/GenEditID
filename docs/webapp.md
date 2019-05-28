# Genome Editing WebApp

## Python Web Frameworks

### Links to articles

- [Wikipedia article on 'Comparison of web frameworks'](https://en.wikipedia.org/wiki/Comparison_of_web_frameworks#Python_2)
- [Django vs Flask vs Pyramid: Choosing a Python Web Framework](https://www.airpair.com/python/posts/django-flask-pyramid)
- [Python Framework Comparison: Django vs. Pyramid](https://www.codementor.io/sheena/django-vs-pyramid-python-framework-comparison-du107yb1c)

### The options

- [Pyramid](https://trypyramid.com/)
- [Django](https://www.djangoproject.com/)
- [Flask](http://flask.pocoo.org/)

### What we use: Pyramid, Chameleon, Deform and Plotly links

- http://docs.pylonsproject.org/projects/pyramid/en/latest/
- http://chameleon.readthedocs.io/en/latest/reference.html
- http://docs.pylonsproject.org/projects/pyramid/en/latest/quick_tutorial/forms.html
- http://deformretail.chromaticleaves.com/
- http://docs.pylonsproject.org/projects/deform/en/latest/api.html
- https://plot.ly/python/reference/


## Our first draft of our WebApp design
![Image of web-design](web-design.jpg)


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
  createdb geneediting
  createuser gene
  psql geneediting
  geneediting=# GRANT CREATE ON SCHEMA public TO "gene";
  geneediting=# GRANT USAGE ON SCHEMA public TO "gene";
  geneediting=# ALTER USER gene WITH PASSWORD 'gene';
  # edit `crispr.yml` to set DATABASE_URI: "postgresql://gene:gene@localhost/geneediting"
  python python/scripts/create_db.py
  ```
- Run your webapp
  ```
  pserve development.ini --reload
  ```
  Go to http://localhost:8080
- Create a new project
- Load project data
  ```
  ./shell/load_project_GEP00001.sh
  ```


## Running the WebApp after first setup

Edit `python/webapp/development.ini` to check you are pointing to the right database

```
cd editID/python/webapp
source wenv/bin/activate
pserve development.ini --reload
```

Go to http://localhost:8080


## Setting up the server for production

https://docs.pylonsproject.org/projects/pyramid/en/latest/tutorials/modwsgi/index.html

Once you have Apache installed, install mod_wsgi.
https://github.com/GrahamDumpleton/mod_wsgi

```bash
ssh bioinf-ge001.cri.camres.org
sudo yum install httpd httpd-devel
sudo yum install python34-devel

sudo su - ge
source venv/bin/activate
pip install mod_wsgi
mod_wsgi-express start-server
```

Verify that the installation worked by pointing your browser at:
http://bioinf-ge001.cri.camres.org:8000/

Create `genome-editing/python/webapp/pyramid.wsgi` with this content
```python
from pyramid.paster import get_app, setup_logging
ini_path = '/home/ge/genome-editing/python/webapp/production.ini'
setup_logging(ini_path)
application = get_app(ini_path, 'main')
```

Configure Apache, with file `vim /etc/httpd/conf.d/genomeediting.conf`
NB. remove specific configuration for the time being otherwise links generated
by Pyramid become absolute and come up as http://localhost:8080/ which is a problem.
Restart Apache
```bash
ssh bioinf-ge001.cri.camres.org
sudo /sbin/service httpd restart
```

Starting the web app:
```
sudo su - ge
source venv/bin/activate
cd genome-editing/python/webapp
mod_wsgi-express start-server pyramid.wsgi --port 8080 --user ge --group users --server-root=/home/ge/
```
http://bioinf-ge001.cri.camres.org:8080/

Restarting the web app:

```bash
ssh bioinf-ge001.cri.camres.org
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

Log file /home/ge/error_log

Go to http://bioinf-ge001.cri.camres.org:8080/

NB. running on port 8080, apache config needs to be fixed.
