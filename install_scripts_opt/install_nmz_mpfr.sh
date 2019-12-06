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
    export LDFLAGS="${LDFLAGS} -L`brew --prefix`/opt/llvm/lib"
fi

## script for the installation of MPFR (needed for Flint)

MPFR_VERSION="4.0.1"
MPFR_URL="http://www.mpfr.org/mpfr-${MPFR_VERSION}/mpfr-${MPFR_VERSION}.tar.gz"
MPFR_SHA256=e650f8723bfc6eca4f222c021db3d5d4cebe2e21c82498329bb9e6815b99c88c

if [ "x$NMZ_PREFIX" != x ]; then
    mkdir -p ${NMZ_PREFIX}
    PREFIX=${NMZ_PREFIX}
else
    PREFIX=${PWD}/local
fi

echo "Installing MPFR..."

mkdir -p ${NMZ_OPT_DIR}/MPFR_source/
cd ${NMZ_OPT_DIR}/MPFR_source
../../download.sh ${MPFR_URL} ${MPFR_SHA256}
if [ ! -d mpfr-${MPFR_VERSION} ]; then
    tar -xvf mpfr-${MPFR_VERSION}.tar.gz
fi
cd mpfr-${MPFR_VERSION}
if [ ! -f config.status ]; then
    ./configure --prefix=${PREFIX} $WITH_GMP
fi
make -j4
make install
