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

## script for the installation of e-antic for the use in libnormaliz
E_ANTIC_VERSION=0.1.2

if [ "x$NMZ_PREFIX" != x ]; then
    mkdir -p ${NMZ_PREFIX}
    PREFIX=${NMZ_PREFIX}
else
    PREFIX=${PWD}/local
fi

if [ "x$NO_OPENMP" != x ]; then
    export BLOCK_OPENMP="--disable-openmp"
fi

echo "Installing E-ANTIC..."

mkdir -p ${NMZ_OPT_DIR}/E-ANTIC_source/
cd ${NMZ_OPT_DIR}/E-ANTIC_source
if [ ! -d e-antic-${E_ANTIC_VERSION} ]; then
    wget http://www.labri.fr/perso/vdelecro/e-antic/e-antic-${E_ANTIC_VERSION}.tar.gz
    tar -xvf e-antic-${E_ANTIC_VERSION}.tar.gz
fi
cd e-antic-${E_ANTIC_VERSION}
if [ ! -f config.status ]; then
    ./configure --prefix=${PREFIX} $WITH_GMP  ${BLOCK_OPENMP} CFLAGS=-I${PREFIX}/include \
              CPPFLAGS="-I${PREFIX}/include -fPIC" \
              LDFLAGS=-L/${PREFIX}/lib
# --enable-flint-devel ## for Flint development version
fi
make -j4
make install

