# Server set up

Information on setting up the server from scratch (actually from a minimal install of Centos 7).

## Postgres

See [postgres.md](the Postgres document).

## Java

Download the JDK from [http://www.oracle.com/technetwork/java/javase/downloads/index.html](Oracle).
Fetch the Linux x64 RPM version.

Install with:

```
sudo yum localinstall -y jdk-8u121-linux-x64.rpm
```

## RPM Packages

There are several things needed to be installed through YUM to build R and its packages. There will also be other tools required for general use.

```
yum install \
    cairo-devel \
    gcc-c++ \
    gcc-gfortran \
    git \
    libcurl-devel \
    openssl-devel \
    pcre-devel \
    xz-devel
```

## R

R can be installed through the Centos package manager (`yum install R`) but it has been decided to install from source.
The installation of R 3.3.3 is made to `/usr/local/R-3.3.3`.

### Base Installation

```
wget https://cran.r-project.org/src/base/R-3/R-3.3.3.tar.gz
tar xfz R-3.3.3.tar.gz

cd R-3.3.3
./configure --prefix=/usr/local/R-3.3.3 --without-x
make
sudo make install
```

### Additional Packages

Some packages are required and to keep them common to this R install, they can be built into R itself.
(This is not possible if we use R from the EPEL RPMs.)

These commands need to be run as root.

```
export PG_INCDIR=/usr/pgsql-9.3/include
export PG_LIBDIR=/usr/pgsql-9.3/lib
/usr/local/R-3.3.3/bin/R
install.packages(c('shiny', 'reshape2', 'ggplot2', 'grofit', 'plotly', 'svglite', 'dplyr', 'RPostgreSQL'),
                 repos="http://mirrors.ebi.ac.uk/CRAN/")
```
