#!/usr/bin/env bash

set -e

echo "::group:e-antic"

source $(dirname "$0")/common.sh

if [ "$GMP_INSTALLDIR" != "" ]; then
    CPPFLAGS="${CPPFLAGS} -I${GMP_INSTALLDIR}/include"
    LDFLAGS="${LDFLAGS} -L${GMP_INSTALLDIR}/lib"
fi

## script for the installation of e-antic for the use in libnormaliz
E_ANTIC_VERSION=0.1.9
E_ANTIC_URL="https://www.labri.fr/perso/vdelecro/e-antic/e-antic-${E_ANTIC_VERSION}.tar.gz"
E_ANTIC_SHA256=f73dad444fd81422d9a1594668615da94405fb6428ce8f2535b58f4cfe88a72d

CONFIGURE_FLAGS="--prefix=${PREFIX}"

if [ "x$NO_OPENMP" != x ]; then
    CONFIGURE_FLAGS="${CONFIGURE_FLAGS} --disable-openmp"
fi

echo "Installing E-ANTIC..."

mkdir -p ${NMZ_OPT_DIR}/E-ANTIC_source/
cd ${NMZ_OPT_DIR}/E-ANTIC_source

../../download.sh ${E_ANTIC_URL} ${E_ANTIC_SHA256}
if [ ! -d e-antic-${E_ANTIC_VERSION} ]; then
    tar -xvf e-antic-${E_ANTIC_VERSION}.tar.gz

fi

cd e-antic-${E_ANTIC_VERSION}

if [ ! -f config.status ]; then
    ./configure ${CONFIGURE_FLAGS}  CFLAGS="${CFLAGS} -I${PREFIX}/include" \
              CPPFLAGS="${CPPFLAGS}" LDFLAGS="${LDFLAGS}"
# --enable-flint-devel ## for Flint development version
fi
make -j4
make install
