# How to...

## Dedicated Host

We have a machine for this project: `bioinf-ge001.cri.camres.org`. This is a virtual machine running Centos 7.

Anne, Chandu, Rich & Ruben have sudo root access to the machine.

Information on setting this up can be found in [Server set up](server_setup.md).

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

#### The options
- [SQLite](https://sqlite.org/) for testing
  - Python3 module to connect to SQLite is called [sqlite3](https://docs.python.org/3.6/library/sqlite3.html#module-sqlite3) and comes by default with Python
- [Postgresql](https://www.postgresql.org/) for production
  - Python3 module to connect to Postgres is called [psycopg2](http://initd.org/psycopg/) and needs to be installed separately (from `requirements.txt`)
- [SQLAlchemy - The Database Toolkit and Object Relational Mapper for Python](http://www.sqlalchemy.org/) gives access to any databases in Python and needs to be installed separately (from `requirements.txt`)
  - [SQLAlchemy 1.1 Documentation](http://docs.sqlalchemy.org/en/rel_1_1/)
- [Alembic](https://bitbucket.org/zzzeek/alembic) is an excellent solution for SQLAlchemy-based systems. It provides a methodical approach and supports auto-generation of migration scripts. See [article](https://www.compose.com/articles/schema-migrations-with-alembic-python-and-postgresql/).

#### Installation instruction
Installation instructions for Postgres on a Centos 7 server can be found [Installing Postgres](postgres.md).


#### Postgres Drivers for R

## R

### Installing R 3.3.2
[https://cran.r-project.org/bin/macosx/](https://cran.r-project.org/bin/macosx/)

R 3.3.2 binary for Mac OS X 10.9 (Mavericks) and higher, signed package. Contains R 3.3.2 framework, R.app GUI 1.68 in 64-bit for Intel Macs, Tcl/Tk 8.6.0 X11 libraries and Texinfo 5.2. The latter two components are optional and can be omitted when choosing "custom install", it is only needed if you want to use the tcltk R package or build package documentation from sources.

### Installing RStudio 1.0.136
[https://www.rstudio.com/products/rstudio/download/](https://www.rstudio.com/products/rstudio/download/)

### R packages dependencies

```R
install.packages(c('shiny', 'reshape2', 'ggplot2', 'grofit', 'plotly', 'svglite', 'dplyr', 'RColorBrewer', 'RSQLite','RPostgreSQL', 'DT'), repos="http://mirrors.ebi.ac.uk/CRAN/")
source("https://bioconductor.org/biocLite.R")
biocLite("ggbio")
```

## Create database schema

Install [dependencies](#dependencies) first.


```bash
source venv/bin/activate
export PYTHONPATH=`pwd`/python
python python/scripts/create_db.py
```

Visualize the SQLite database using [DbVisualizer](http://www.dbvis.com/).

### Create database schema on dedicated server

- Edit configuration file `python/dnascissors/crispr.yml` file and use `DATABASE_URI: "postgresql://gene:gene@bioinf-ge001.cri.camres.org/geneediting"`
- Run `python/scripts/create_db.py` script to create DB schema

## Load data

Four scripts for loading into the tracking database:
- project data and plate layouts,
- protein abundance (ICW channels),
- cell growths (Incucyte) and
- variant results from NGS analysis.

```bash
source venv/bin/activate
export PYTHONPATH=`pwd`/python
python python/scripts/load_layout.py --layout=data/20170127_GEP00001/20170118_GEP00001.xlsx
python python/scripts/load_protein_abundance.py --plateid=GEP00001_01 --file=data/20170127_GEP00001/GEP00001_01_ICW.csv
python python/scripts/load_cell_growth.py --plateid=GEP00001_01 --file=data/20170127_GEP00001/GEP00001_01_incu.txt
python python/scripts/load_variant_results.py --file=data/20170127_GEP00001/GEP00001_NGS_IndelsResults.csv
```

One script to load all files associated with project GEP00001

```bash
shell/load_project_GEP00001.sh
```

### Load data into database on dedicated server

- Edit configuration file `python/dnascissors/crispr.yml` file and use `DATABASE_URI: "postgresql://gene:gene@bioinf-ge001.cri.camres.org/geneediting"`
- Run `shell/load_project_GEP00001.sh` script to load all data associated to GEP00001 project

### Database migration

Using alembic http://alembic.zzzcomputing.com/en/latest/index.html

* Creating an environment
  ```bash
  source venv/bin/activate
  pip install alambic
  cd python
  alembic init dbmigration
  ```
* Creating and running a migration script automatically
  ```bash
  # modify model.py and automatically generate the changes
  alembic revision --autogenerate -m 'add columns to project'

  alembic upgrade head
  ```


## NGS analysis

See these files for more information:
- [NGS Pipeline](ngs-pipeline.md)
- [Run the Amplicon Sequencing Pipeline](ngs-run-pipeline.md)
- [How to visualise NGS data and variants?](ngs-data-vis.md)
- [NGS Downstream processing and plotting](ngs-downstream.md)

## Results in Shiny App

### R script for plotting from DB
- First, install R packages (see list above) in RStudio
- run script `r/scripts/genome_editing.r` to plot protein abundance and clone growth curve.

### How to build a shiny app
- [http://shiny.rstudio.com/tutorial/](http://shiny.rstudio.com/tutorial/)
- [http://shiny.rstudio.com/tutorial/lesson1/](http://shiny.rstudio.com/tutorial/lesson1/)

### How to run the shiny app in RStudio
- Installation instruction
  - in RStudio, File > New Project and select directory of the git repo
  - install these R packages in RStudio (see above for the full list)
  - Run the app
  ```R
  shiny::runApp('r/shinyapp', port=4700, host='0.0.0.0')
  ```

### How to run the shiny app from R or on the server.

```bash
R -e "shiny::runApp('r/shinyapp', port=4700, host='0.0.0.0')"
```

## Loaders and Results in WebApp

Our first draft of our WebApp design
![Image of web-design](web-design.jpg)

For more information see [CRISPR WebApp](webapp.md).
