#!/usr/bin/env bash

set -e

WITH_GMP=""
if [ "$GMP_INSTALLDIR" != "" ]; then
  WITH_GMP="--with-gmp=$GMP_INSTALLDIR"
fi

if [ "x$NMZ_OPT_DIR" = x ]; then 
    export NMZ_OPT_DIR=${PWD}/nmz_opt_lib
    mkdir -p ${NMZ_OPT_DIR}
fi

if [ "x$NMZ_COMPILER" != x ]; then
    export CXX=$NMZ_COMPILER
elif [[ $OSTYPE == darwin* ]]; then
    export CXX=clang++
    export PATH="`brew --prefix`/opt/llvm/bin/:$PATH"
    export LDFLAGS="-L`brew --prefix`/opt/llvm/lib"
fi

## script for the installation of Flint for the use in libnormaliz
## including the installation of MPFR (needed for Flint)

FLINT_VERSION="2.5.2"
MPFR_VERSION="4.0.1"

if [ "x$NMZ_PREFIX" != x ]; then
    mkdir -p ${NMZ_PREFIX}
    PREFIX=${NMZ_PREFIX}
else
    PREFIX=${PWD}/local
fi

echo "Installing MPFR..."

mkdir -p ${NMZ_OPT_DIR}/MPFR_source/
cd ${NMZ_OPT_DIR}/MPFR_source
if [ ! -d mpfr-${MPFR_VERSION} ]; then
    wget http://www.mpfr.org/mpfr-${MPFR_VERSION}/mpfr-${MPFR_VERSION}.tar.gz
    tar -xvf mpfr-${MPFR_VERSION}.tar.gz
fi
cd mpfr-${MPFR_VERSION}
./configure --prefix=${PREFIX} $WITH_GMP
make -j4
make install

echo "Installing FLINT..."

mkdir -p ${NMZ_OPT_DIR}/Flint_source/
cd ${NMZ_OPT_DIR}/Flint_source
if [ ! -d flint-${FLINT_VERSION} ]; then
    wget http://www.flintlib.org/flint-${FLINT_VERSION}.tar.gz
    tar -xvf flint-${FLINT_VERSION}.tar.gz
fi
cd flint-${FLINT_VERSION}
if [ ! -f Makefile ]; then
    ./configure --prefix=${PREFIX} --with-mpfr=${PREFIX} $WITH_GMP $EXTRA_FLINT_FLAGS
fi
# patch to avoid PIE clash in Ubuntu >= 16-10
if [[ $OSTYPE == "linux-gnu" ]]; then
    sed -i s/"-Wl,"// Makefile.subdirs
fi
make -j4 verbose
make install
