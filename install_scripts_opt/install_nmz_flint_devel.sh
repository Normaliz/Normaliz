#!/usr/bin/env bash

set -e

WITH_GMP=""
if [ "$GMP_INSTALLDIR" != "" ]; then
  WITH_GMP="--with-gmp=$GMP_INSTALLDIR"
fi

source $(dirname "$0")/common.sh

## script for the installation of Flint for the use in libnormaliz

FLINT_BRANCH=trunk
FLINT_COMMIT=5be51316aff24f6b04edcfe3dd4388bdff2934db

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

if [[  $CXX == clang* ]]; then        
    ## To prevent a syntax error in clang with -fopenmp
    sed -i.orig 's/#pragma omp/\/\/ #pragma omp/;' $PREFIX/include/flint/ulong_extras.h
fi
