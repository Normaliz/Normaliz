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

FLINT_BRANCH=trunk
FLINT_COMMIT=a4a87533fc113875bb0f1ed0156ea1319de3051c
MPFR_VERSION="4.0.0"

#PREFIX=${NMZ_OPT_DIR}
PREFIX=${PWD}/local


echo "Installing MPFR..."

mkdir -p ${NMZ_OPT_DIR}/MPFR_source/
cd ${NMZ_OPT_DIR}/MPFR_source
if [ ! -d mpfr-${MPFR_VERSION} ]; then
    curl -O http://www.mpfr.org/mpfr-${MPFR_VERSION}/mpfr-${MPFR_VERSION}.tar.gz
    tar -xvf mpfr-${MPFR_VERSION}.tar.gz
fi
cd mpfr-${MPFR_VERSION}
if [ ! -f config.status ]; then
    ./configure --prefix=${PREFIX} $WITH_GMP
fi
make -j4
make install

echo "Installing FLINT..."

mkdir -p ${NMZ_OPT_DIR}/Flint_source/
cd ${NMZ_OPT_DIR}/Flint_source
if [ ! -d flint2 ]; then
    git clone --branch=${FLINT_BRANCH} --single-branch https://github.com/wbhart/flint2.git
    (cd flint2 && git checkout ${FLINT_COMMIT})
fi
cd flint2
if [ ! -f Makefile ]; then
    ./configure --prefix=${PREFIX} --with-mpfr=${PREFIX} $WITH_GMP $EXTRA_FLINT_FLAGS
fi
make -j4 verbose
make install

if [[  $OSTYPE == darwin* ]]; then        
    ## To prevent a syntax error in clang with -fopenmp
    sed -i.orig 's/#pragma omp/\/\/ #pragma omp/;' local/include/flint/ulong_extras.h
fi
