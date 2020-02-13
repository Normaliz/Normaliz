#!/usr/bin/env bash

set -e

WITH_GMP=""
if [ "$GMP_INSTALLDIR" != "" ]; then
  WITH_GMP="--with-gmp=$GMP_INSTALLDIR"
fi

source $(dirname "$0")/common.sh

## script for the installation of e-antic for the use in libnormaliz
E_ANTIC_VERSION=0.1.3b0

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

