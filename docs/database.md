# GenEditID database


## Installing packages on a linux server

### on RHEL 7
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

### On ubuntu
(As root)
```
sudo apt install libpq-dev postgresql-client
  # the libpq-dev is required to install the R package RPostgreSQL, which otherwise throws an error:
  # RS-PostgreSQL.h:23:26: fatal error: libpq-fe.h: No such file or directory compilation terminated.
  # /usr/lib/R/etc/Makeconf:134: recipe for target 'RS-PQescape.o' failed
sudo apt install postgresql postgresql-contrib
```

## Configuring access to the database

(As postgres)
```
sudo -i -u postgres
```

Edit `/var/lib/pgsql/9.3/data/pg_hba.conf` to add access to the (to be created) database `geneditid`. Add the following:
```
# TYPE  DATABASE        USER        ADDRESS        METHOD
local    geneditid    postgres                   ident
local    geneditid    gene                       md5
host     geneditid    gene        127.0.0.1/8    md5
host     geneditid    gene        ::1/128        md5
```

Edit `/var/lib/pgsql/9.3/data/postgresql.conf` to change the `listen_addresses` to "`*`" (line 59).

(As root)
```
systemctl enable postgresql-9.3
systemctl start postgresql-9.3
```

(As postgres again)
```
createdb geneditid
createuser gene
psql geneditid
```

(This is now in the Postgres client.)
```
GRANT CREATE ON SCHEMA public TO "gene";
GRANT USAGE ON SCHEMA public TO "gene";
ALTER USER gene WITH PASSWORD 'gene';
```

## Accessing the database from elsewhere

```
psql -h <HOST_URL> -U gene geneditid
```

The YAML setting for accessing this database in the code is:
```
DATABASE_URI: "postgresql://gene:gene@<HOST_URL>/geneditid"
```

## Python dependencies
- [Postgresql](https://www.postgresql.org/) for production
  - Python3 module to connect to Postgres is called [psycopg2](http://initd.org/psycopg/) and needs to be installed separately (from `requirements.txt`)
- [SQLAlchemy - The Database Toolkit and Object Relational Mapper for Python](http://www.sqlalchemy.org/) gives access to any databases in Python and needs to be installed separately (from `requirements.txt`)
  - [SQLAlchemy 1.1 Documentation](http://docs.sqlalchemy.org/en/rel_1_1/)
- [Alembic](https://bitbucket.org/zzzeek/alembic) is an excellent solution for SQLAlchemy-based systems. It provides a methodical approach and supports auto-generation of migration scripts. See [article](https://www.compose.com/articles/schema-migrations-with-alembic-python-and-postgresql/).


## Create database schema

```
git clone https://github.com/GenEditID/editID.git
cd editID/
python3 -m venv venv
source venv/bin/activate
pip install -r python/requirements.txt

python python/scripts/create_db.py
```

Access the database using [DbVisualizer](http://www.dbvis.com/).

View the current [database schema](db_diagram.pdf).


## Create database schema on dedicated server

- Edit configuration file `python/dnascissors/crispr.yml` file and use `DATABASE_URI: "postgresql://gene:gene@<HOST_URL>/geneditid"`
- Run `python/scripts/create_db.py` script to create DB schema
