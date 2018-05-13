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

## script for the installation of ARB for the use in libnormaliz

ARB_VERSION="2.13.0"

#PREFIX=${NMZ_OPT_DIR}
PREFIX=${PWD}/local

echo "Installing ARB..."

mkdir -p ${NMZ_OPT_DIR}/ARB_source/
cd ${NMZ_OPT_DIR}/ARB_source
if [ ! -d arb-${ARB_VERSION} ]; then
    curl -L -o arb-${ARB_VERSION}.tar.gz -O https://github.com/fredrik-johansson/arb/archive/${ARB_VERSION}.tar.gz
    tar -xvf arb-${ARB_VERSION}.tar.gz
fi
cd arb-${ARB_VERSION}
# (In particular on Mac OS X, make sure that our version of MPFR comes
# first in the -L search path, not the one from LLVM or elsewhere.
# ARB's configure puts it last.)
## export LDFLAGS="-L${NMZ_OPT_DIR}/lib ${LDFLAGS}"
export LDFLAGS="-L${PREFIX}/lib ${LDFLAGS}"
if [ ! -f Makefile ]; then
    ./configure --prefix=${PREFIX} $WITH_GMP --with-flint="${PREFIX}" \
                --with-mpfr="${PREFIX}"
fi
make -j4 verbose
make install
