# Genome Editing project

[![Gitter](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/genome-editing/Lobby?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

## Technologies used: Postgres, Python and R

* DB schema to store all data
* Database: Postgresql / SQLite
* DB schema/mapper: DDLutil [https://db.apache.org/ddlutils/](https://db.apache.org/ddlutils/) (XML to db) or SQLAlchemy [https://www.sqlalchemy.org/](https://www.sqlalchemy.org/) (Object Relational Mapper)
* Python [https://www.python.org/](https://www.python.org/)
* R [https://www.r-project.org/](https://www.r-project.org/) and shiny [https://www.rstudio.com/products/shiny/](https://www.rstudio.com/products/shiny/) (Easy web applications in R)
* Editors:
  - PyCharm [https://www.jetbrains.com/pycharm/](https://www.jetbrains.com/pycharm/)
  - Eclipse [http://www.eclipse.org/](http://www.eclipse.org/) with Python plugin PyDev [http://www.pydev.org/](http://www.pydev.org/)
  - RStudio [https://www.rstudio.com/](https://www.rstudio.com/)
  - Atom [https://atom.io/](https://atom.io/) or any others you like!

## Project description

See [docs/project_description.md](docs/project_description.md)

* three modules: growth, western and NGS
* three types of files so far:
  - PlatesLayout_270117.csv
  - 240117_ICW_Plate1.csv (possible to have 8 colours)
  - 180117_MCF7_luc-straw_clone3_CRISPR_plate1.txt
* one shiny app in `r/shinyapp` directory with instruction on how to run it in [docs/howto.md](docs/howto.md#shiny-app)
* database with instruction on how to create it in [docs/howto.md](docs/howto.md#create-database)

## People involved

* Ruben Alvarez [@rubenalv](https://github.com/rubenalv)
* Rich Bowers [@rich7409](https://github.com/rich7409)
* Chandu Chilamakuri [@chilamakuricsreddy](https://github.com/chilamakuricsreddy)
* Anne Pajon [@pajanne](https://github.com/pajanne)

## Planning discussion and TODO list

See [planning](planning.md)

## How to...

See [How to...](docs/howto.md) install and run this project

## Need help with Git?

See [Getting started with Git](docs/git-help.md)
