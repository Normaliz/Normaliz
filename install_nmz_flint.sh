#!/bin/sh -e
if [ "x$NMZ_OPT_DIR" = x ]; then 
    export NMZ_OPT_DIR=${PWD}/nmz_opt_lib
    mkdir -p ${NMZ_OPT_DIR}
fi

## script for the installation of Flint for the use in libnormaliz
## including the installation of MPFR (needed for Flint)

FLINT_VERSION="2.5.2"
MPFR_VERSION="4.0.0"

PREFIX=${NMZ_OPT_DIR}

mkdir -p ${NMZ_OPT_DIR}/MPFR_source/
cd ${NMZ_OPT_DIR}/MPFR_source
wget -N http://www.mpfr.org/mpfr-current/mpfr-${MPFR_VERSION}.tar.gz
tar -xvf mpfr-${MPFR_VERSION}.tar.gz
cd mpfr-${MPFR_VERSION}
./configure --prefix=${PREFIX}
make -j4
make install

mkdir -p ${NMZ_OPT_DIR}/Flint_source/
cd ${NMZ_OPT_DIR}/Flint_source
wget -4 -N http://www.flintlib.org/flint-${FLINT_VERSION}.tar.gz
tar -xvf flint-${FLINT_VERSION}.tar.gz
cd flint-${FLINT_VERSION}
./configure --prefix=${PREFIX} --with-mpfr=${PREFIX}
make -j4
make install