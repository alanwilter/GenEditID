# Genome Editing Database

Notes on installing Postgres on RHEL 7.

For this, I have installed Postgres 9.3 simply because the Clarity servers are running 9.3 and it makes Anne and my life easier to have the same version. Postgres is on 9.6 but I very much doubt we'll be using features in the later one that aren't in the earlier.

The shared instance is running on bioinf-srv003 until we get a dedicated VM.


## Installing Packages

(As root)
```
yum localinstall https://download.postgresql.org/pub/repos/yum/9.3/redhat/rhel-7-x86_64/pgdg-centos93-9.3-3.noarch.rpm
yum install postgresql93 postgresql93-server postgresql93-devel

passwd postgres
```

(Set the password to something usable. It won't be used much.)

```
su - postgres
/usr/pgsql-9.3/bin/initdb
exit
```


## On Ubuntu
(As root)
```
sudo apt install libpq-dev postgresql-client
  # the libpq-dev is required to install the R package RPostgreSQL, which otherwise throws an error:
  # RS-PostgreSQL.h:23:26: fatal error: libpq-fe.h: No such file or directory compilation terminated.
  # /usr/lib/R/etc/Makeconf:134: recipe for target 'RS-PQescape.o' failed
sudo apt install postgresql postgresql-contrib
```


## Configuring Access to the Database

(As postgres)

```
sudo -i -u postgres
```

Edit `/var/lib/pgsql/9.3/data/pg_hba.conf` to add access to the (to be created) database `geneediting`. Add the following:

```
# TYPE  DATABASE        USER        ADDRESS        METHOD
local    geneediting    postgres                   ident
local    geneediting    gene                       md5
host     geneediting    gene        127.0.0.1/8    md5
host     geneediting    gene        ::1/128        md5
host     geneediting    gene        10.20.0.0/16   md5
```

Edit `/var/lib/pgsql/9.3/data/postgresql.conf` to change the `listen_addresses` to "`*`" (line 59).

(As root)

```
systemctl enable postgresql-9.3
systemctl start postgresql-9.3
```

(As postgres again)

```
createdb geneediting
createuser gene
psql geneediting
```

(This is now in the Postgres client.)

```
GRANT CREATE ON SCHEMA public TO "gene";
GRANT USAGE ON SCHEMA public TO "gene";
ALTER USER gene WITH PASSWORD 'gene';
```


## Accessing the Database from elsewhere

```
psql -h bioinf-ge001.cri.camres.org -U gene geneediting
```

The YAML setting for accessing this database in the code is:

```
DATABASE_URI: "postgresql://gene:gene@bioinf-ge001.cri.camres.org/geneediting"
```


## Python dependencies
- [Postgresql](https://www.postgresql.org/) for production
  - Python3 module to connect to Postgres is called [psycopg2](http://initd.org/psycopg/) and needs to be installed separately (from `requirements.txt`)
- [SQLAlchemy - The Database Toolkit and Object Relational Mapper for Python](http://www.sqlalchemy.org/) gives access to any databases in Python and needs to be installed separately (from `requirements.txt`)
  - [SQLAlchemy 1.1 Documentation](http://docs.sqlalchemy.org/en/rel_1_1/)
- [Alembic](https://bitbucket.org/zzzeek/alembic) is an excellent solution for SQLAlchemy-based systems. It provides a methodical approach and supports auto-generation of migration scripts. See [article](https://www.compose.com/articles/schema-migrations-with-alembic-python-and-postgresql/).


## Create database schema

```shell
source venv/bin/activate
export PYTHONPATH=`pwd`/python
python python/scripts/create_db.py
```

Access the database using [DbVisualizer](http://www.dbvis.com/).

View the [database schema](db_diagram.pdf).


## Create database schema on dedicated server

- Edit configuration file `python/dnascissors/crispr.yml` file and use `DATABASE_URI: "postgresql://gene:gene@bioinf-ge001.cri.camres.org/geneediting"`
- Run `python/scripts/create_db.py` script to create DB schema


## Database migration

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
