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

## script for the installation of ANTIC for the use in libnormaliz

ANTIC_BRANCH=trunk
ANTIC_COMMIT=36f19eabcd7c2051fe3ed9b5ff54ba5816d7aba7
PREFIX=${NMZ_OPT_DIR}

echo "Installing ANTIC..."

mkdir -p ${NMZ_OPT_DIR}/ANTIC_source/
cd ${NMZ_OPT_DIR}/ANTIC_source
if [ ! -d antic ]; then
    git clone --branch=${ANTIC_BRANCH} --single-branch https://github.com/wbhart/antic.git
    (cd antic && git checkout ${ANTIC_COMMIT})
fi
cd antic
# (In particular on Mac OS X, make sure that our version of MPFR comes
# first in the -L search path, not the one from LLVM or elsewhere.
# ANTIC's configure puts it last.)
export LDFLAGS="-L${NMZ_OPT_DIR}/lib ${LDFLAGS}"
if [ ! -f Makefile ]; then
    ./configure --prefix=${PREFIX} $WITH_GMP --with-flint="${NMZ_OPT_DIR}" \
                --with-mpfr="${NMZ_OPT_DIR}"
fi
make -j4 verbose
make install
