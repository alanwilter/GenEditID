# How to...

## Dedicated Host

We have a machine for this project: `bioinf-ge001.cri.camres.org`. This is a virtual machine running Centos 7.

Anne, Chandu, Rich & Ruben have sudo root access to the machine.

Information on setting this up can be found in [Server set up](server_setup.md).


## Python

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


## Database

### Dependencies
- [Postgresql](https://www.postgresql.org/) for production
  - Python3 module to connect to Postgres is called [psycopg2](http://initd.org/psycopg/) and needs to be installed separately (from `requirements.txt`)
- [SQLAlchemy - The Database Toolkit and Object Relational Mapper for Python](http://www.sqlalchemy.org/) gives access to any databases in Python and needs to be installed separately (from `requirements.txt`)
  - [SQLAlchemy 1.1 Documentation](http://docs.sqlalchemy.org/en/rel_1_1/)
- [Alembic](https://bitbucket.org/zzzeek/alembic) is an excellent solution for SQLAlchemy-based systems. It provides a methodical approach and supports auto-generation of migration scripts. See [article](https://www.compose.com/articles/schema-migrations-with-alembic-python-and-postgresql/).

### Installation instruction
Installation instructions for Postgres on a Centos 7 server can be found [Installing Postgres](database.md).

### Create database schema

```shell
source venv/bin/activate
export PYTHONPATH=`pwd`/python
python python/scripts/create_db.py
```

Access the database using [DbVisualizer](http://www.dbvis.com/).

View the [database schema](db_diagram.pdf).

### Create database schema on dedicated server

- Edit configuration file `python/dnascissors/crispr.yml` file and use `DATABASE_URI: "postgresql://gene:gene@bioinf-ge001.cri.camres.org/geneediting"`
- Run `python/scripts/create_db.py` script to create DB schema


## Load data

Six scripts for loading data into the database:
- reference data called `python/scripts/load_ref_data.py`,
- project data and plate layouts called `python/scripts/load_layout.py`,
- protein abundance (ICW channels) called `python/scripts/load_protein_abundance.py`,
- cell growths (Incucyte) called `python/scripts/load_cell_growth.py`,
- variant results from NGS analysis called `python/scripts/load_variant_results.py` and
- mutation summary called `python/scripts/load_mutation_summary.py`

One script to load all files associated with project GEP00001:
```
shell/load_project_GEP00001.sh
```

### Load data into database on dedicated server

- Edit configuration file `python/dnascissors/crispr.yml` file and use `DATABASE_URI: "postgresql://gene:gene@bioinf-ge001.cri.camres.org/geneediting"`
- Run `shell/load_project_GEP00001.sh` script to load all data associated to GEP00001 project


## NGS analysis

See these files for more information:
- [NGS Pipeline](ngs-pipeline.md)
- [Run the Amplicon Sequencing Pipeline](ngs-run-pipeline.md)
- [How to visualise NGS data and variants?](ngs-data-vis.md)
- [NGS Downstream processing and plotting](ngs-downstream.md)


## Loaders and Results in WebApp

For more information see [WebApp](webapp.md).


## R (legacy)

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

### Results in legacy Shiny App

#### R script for plotting from DB
- First, install R packages (see list above) in RStudio
- run script `r/scripts/genome_editing.r` to plot protein abundance and clone growth curve.

#### How to build a shiny app
- [http://shiny.rstudio.com/tutorial/](http://shiny.rstudio.com/tutorial/)
- [http://shiny.rstudio.com/tutorial/lesson1/](http://shiny.rstudio.com/tutorial/lesson1/)

#### How to run the shiny app in RStudio
- Installation instruction
  - in RStudio, File > New Project and select directory of the git repo
  - install these R packages in RStudio (see above for the full list)
  - Run the app
  ```R
  shiny::runApp('r/shinyapp', port=4700, host='0.0.0.0')
  ```

#### How to run the shiny app from R or on the server.

```bash
R -e "shiny::runApp('r/shinyapp', port=4700, host='0.0.0.0')"
```
