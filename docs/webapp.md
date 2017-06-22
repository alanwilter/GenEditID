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


## getting started

See file in `python/webapp/README.txt`

- Change directory into your newly created project.
```bash
    cd webapp
```
- Create a Python virtual environment.
```bash
    python3 -m venv env
```
- Upgrade packaging tools.
```bash
    env/bin/pip install --upgrade pip setuptools
```
- Install the project in editable mode with its testing requirements.
```bash
    env/bin/pip install -e ".[testing]"
```
- Configure the database.
```bash
    env/bin/initialize_webapp_db development.ini
```
- Run your project's tests.
```bash
    env/bin/pytest
```
- Run your project.
```bash
    env/bin/pserve development.ini --reload
```

Go to http://localhost:8080
