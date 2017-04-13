# How to...

## Shiny App

### Installing R 3.3.2
[https://cran.r-project.org/bin/macosx/](https://cran.r-project.org/bin/macosx/)

R 3.3.2 binary for Mac OS X 10.9 (Mavericks) and higher, signed package. Contains R 3.3.2 framework, R.app GUI 1.68 in 64-bit for Intel Macs, Tcl/Tk 8.6.0 X11 libraries and Texinfo 5.2. The latter two components are optional and can be ommitted when choosing "custom install", it is only needed if you want to use the tcltk R package or build package documentation from sources.

### Installing RStudio 1.0.136
[https://www.rstudio.com/products/rstudio/download/](https://www.rstudio.com/products/rstudio/download/)

### How to build a shiny app
- [http://shiny.rstudio.com/tutorial/](http://shiny.rstudio.com/tutorial/)
- [http://shiny.rstudio.com/tutorial/lesson1/](http://shiny.rstudio.com/tutorial/lesson1/)

### How to run the shiny app
- Installation instruction
  - in RStudio, File > New Project and select directory of the git repo
  - install these packages in RStudio
  ```R
  install.packages(c('shiny', 'reshape2', 'ggplot2', 'grofit', 'plotly'))
  ```
- Run the app
```R
source('r/shinyapp/global.R')
runApp('r/shinyapp/.')
```
- Choose CSV File, click on 'Browse...' and select all 13 data files in `data/20170127_Experiment0001/` directory.

## Dependencies

### Python
- Install Python3, check the version:
```
python3 --version
Python 3.6
```
- Get help from the documentation [https://docs.python.org/3.6/](https://docs.python.org/3.6/)
- Create a virtual environment to install Python libraries needed:
```bash
python3 -m venv venv
# activate your virtual environment
source venv/bin/activate
# install all requirements
pip install -r python/requirements.txt
```

### Database

- [SQLite](https://sqlite.org/) for testing
  - Python3 module to connect to SQLite is called [sqlite3](https://docs.python.org/3.6/library/sqlite3.html#module-sqlite3) and comes by default with Python
- [Postgresql](https://www.postgresql.org/) for production
  - Python3 module to connect to Postgres is called [psycopg2](http://initd.org/psycopg/) and needs to be installed separately (from `requirements.txt`)
- [SQLAlchemy - The Database Toolkit and Object Relational Mapper for Python](http://www.sqlalchemy.org/) gives access to any databases in Python and needs to be installed separately (from `requirements.txt`)
  - [SQLAlchemy 1.1 Documentation](http://docs.sqlalchemy.org/en/rel_1_1/)
- [Alembic](https://bitbucket.org/zzzeek/alembic) is an excellent solution for SQLAlchemy-based systems. It provides a methodical approach and supports auto-generation of migration scripts. See [article](https://www.compose.com/articles/schema-migrations-with-alembic-python-and-postgresql/).

Installation instructions for Postgres on a Centos 7 server can be found [in this separate document](postgres.md).


#### Postgres Drivers for R

For R, you will need to install the development files for Postgres:

  yum install postgresql-9.3-devel

Then before starting R, you need to set variables as to where the development files are:

```
export PG_INCDIR=/usr/pgsql-9.3/include
export PG_LIBDIR=/usr/pgsql-9.3/lib
```

Then start R and install:

```
install.packages(c('RPostgreSQL'))
```

## Create database tables

Install [dependencies](#dependencies) first.


```bash
source venv/bin/activate
export PYTHONPATH=`pwd`/python
python python/scripts/create_db.py
```

Visualize the SQLite database using [DbVisualizer](http://www.dbvis.com/).

## Load data

Three scripts for loading the plate layouts, ICW channels and Incucyte growth tracking.

```bash
source venv/bin/activate
export PYTHONPATH=`pwd`/python
python python/scripts/loadLayoutCSV.py
python python/scripts/loadICWCSV.py
python python/scripts/loadIncucyteCSV.py
```

### R script for plotting from DB
- install these packages in RStudio
```R
install.packages(c('dplyr', 'RSQLite', 'ggplot2', 'grofit'))
```
- run script `r/scripts/genome_editing.r` to plot protein abundance and clone growth curve.
