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
source('shinyapp/global.R')
runApp('shinyapp/.')
```
- Choose CSV File, click on 'Browse...' and select all 13 data files in `data` directory.

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
pip install -r requirements.txt
```

### Database

- [SQLite](https://sqlite.org/) for testing
  - Python3 module to connect to SQLite is called [sqlite3](https://docs.python.org/3.6/library/sqlite3.html#module-sqlite3) and comes by default with Python
- [Postgresql](https://www.postgresql.org/) for production
  - Python3 module to connect to Postgres is called [psycopg2](http://initd.org/psycopg/) and needs to be installed separately (from `requirements.txt`)
- [SQLAlchemy - The Database Toolkit and Object Relational Mapper for Python](http://www.sqlalchemy.org/) gives access to any databases in Python and needs to be installed separately (from `requirements.txt`)
  - [SQLAlchemy 1.1 Documentation](http://docs.sqlalchemy.org/en/rel_1_1/)
- [Alembic](https://bitbucket.org/zzzeek/alembic) is an excellent solution for SQLAlchemy-based systems. It provides a methodical approach and supports auto-generation of migration scripts. See [article](https://www.compose.com/articles/schema-migrations-with-alembic-python-and-postgresql/).

## Create database

Install [dependencies](#dependencies) first.


```bash
source venv/bin/activate
python create_db.py
```

Visualize the SQLite database using [DbVisualizer](http://www.dbvis.com/).

## Load data

Next step... load `data/PlatesLayout_270117.csv` into the database.


### R script for plotting from DB
- install these packages in RStudio
```R
install.packages(c('dplyr', 'RSQLite'))
```
