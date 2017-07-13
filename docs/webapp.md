# CRISPR WebApp

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


## Getting started

- Change directory into your newly created project.
```bash
    cd python/webapp
```
- Create a Python virtual environment.
```bash
    python3 -m venv wenv
```
- Activate the virtual environment for the webapp
```bash
    source wenv/bin/activate
```
- Upgrade packaging tools.
```bash
    pip install --upgrade pip setuptools
```
- Install dependencies from the top `python/` directory
```bash
    cd ..
    pip install -r requirements.txt
```
- Setting up a local postgres database
```bash
    createdb geneediting
    createuser gene
    psql geneediting
    geneediting=# GRANT CREATE ON SCHEMA public TO "gene";
    geneediting=# GRANT USAGE ON SCHEMA public TO "gene";
    geneediting=# ALTER USER gene WITH PASSWORD 'gene';
    # edit `crispr.yml` to set DATABASE_URI: "postgresql://gene:gene@localhost/geneediting"
    python python/scripts/create_db.py
    ./shell/load_project_GEP00001.sh
    ./shell/load_project_GEP00002.sh
    ./shell/load_project_GEP00003.sh
```
- Run your project.
```bash
    cd webapp
    pserve development.ini --reload
```

Go to http://localhost:8080/project

## Running the webapp after first setup

Edit `python/webapp/development.ini` to check you are pointing to the right database

```bash
cd python/webapp
source wenv/bin/activate
pserve development.ini --reload
```

Go to http://localhost:8080/project
