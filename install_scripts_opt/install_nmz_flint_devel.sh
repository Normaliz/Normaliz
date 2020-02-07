#!/usr/bin/env bash

set -e

source $(dirname "$0")/common.sh

CONFIGURE_FLAGS="--prefix=${PREFIX}"
if [ "$GMP_INSTALLDIR" != "" ]; then
    CONFIGURE_FLAGS="${CONFIGURE_FLAGS} --with-gmp=$GMP_INSTALLDIR"
fi

# unfortunately flint does not support standard env vars like CPPFLAGS, LDFLAGS;
# but luckily it has some vars of its own that we can use for the same purpose
export EXTRA_INC_DIRS="${PREFIX}/include"
export EXTRA_LIB_DIRS="${PREFIX}/lib"

## script for the installation of Flint for the use in libnormaliz

FLINT_BRANCH=trunk
FLINT_COMMIT=5be51316aff24f6b04edcfe3dd4388bdff2934db

echo "Installing FLINT..."

# download & extract
mkdir -p ${NMZ_OPT_DIR}/Flint_source/
cd ${NMZ_OPT_DIR}/Flint_source
if [ ! -d flint2 ]; then
    git clone --branch=${FLINT_BRANCH} --single-branch https://github.com/wbhart/flint2.git
    (cd flint2 && git checkout ${FLINT_COMMIT})
fi

# configure & compile
cd flint2
if [ ! -f Makefile ]; then
    ./configure ${CONFIGURE_FLAGS} $EXTRA_FLINT_FLAGS
fi
make -j4 verbose
make install

if [[  $CXX == clang* ]]; then        
    ## To prevent a syntax error in clang with -fopenmp
    sed -i.orig 's/#pragma omp/\/\/ #pragma omp/;' $PREFIX/include/flint/ulong_extras.h
fi
