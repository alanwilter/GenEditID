# Server set up

Information on setting up the server from scratch (actually from a minimal install of Centos 7).

## Postgres

See [database.md](the Postgres document).

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
    libxml2-devel \
    openssl-devel \
    pcre-devel \
    screen \
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

These commands need to be run as root. The last step ("ggbio") takes a long time with all its dependencies.

```
export PG_INCDIR=/usr/pgsql-9.3/include
export PG_LIBDIR=/usr/pgsql-9.3/lib

/usr/local/R-3.3.3/bin/R

install.packages(c('shiny', 'reshape2', 'ggplot2', 'grofit', 'plotly', 'svglite', 'dplyr', 'RPostgreSQL', 'DT', 'XML'),
                 repos="http://mirrors.ebi.ac.uk/CRAN/")

source("http://bioconductor.org/biocLite.R")
biocLite("ggbio")
```

# Data Directory

The `bioinf-ge001` server is a small VM with not a lot of storage. Its data directory has therefore
been put onto `bioinf-nfs001` and NFS mounted (with _autofs_).

#### bioinf-nfs001

Create a user matching the "ge" user on bioinf-ge001:

```
useradd -c "Genome Editing" -d /data/ge -m -u 1002 -g sec_bioinf-nfs001-data_rw -G users -N  ge
chmod 775 /data/ge
passwd ge
```

Add the `/data/ge` directory to the NFS exports by adding this to `/etc/exports`:

```
/data/ge     bioinf-ge001(rw,fsid=0,no_subtree_check)
```

Can then export without restart with:

```
exportfs -ra
```

#### bioinf-ge001

Set up _autofs_ to mount `/data/ge` when required over NFS:

```
yum install autofs
```

Then in `/etc/auto.master.d` create two files:

##### nfs.autofs

```
/data /etc/auto.master.d/autofs.nfs
```

##### autofs.nfs

```
ge -fstype=nfs,nosuid bioinf-nfs001.cri.camres.org:/data/ge

```
