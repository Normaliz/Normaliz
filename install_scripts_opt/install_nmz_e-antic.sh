#!/usr/bin/env bash

set -e

source $(dirname "$0")/common.sh

if [ "$GMP_INSTALLDIR" != "" ]; then
    CPPFLAGS="${CPPFLAGS} -I${GMP_INSTALLDIR}/include"
    LDFLAGS="${LDFLAGS} -L${GMP_INSTALLDIR}/lib"
fi

## script for the installation of e-antic for the use in libnormaliz
E_ANTIC_VERSION=1.0.0-rc.11
E_ANTIC_URL="https://github.com/flatsurf/e-antic/releases/download/${E_ANTIC_VERSION}/e-antic-${E_ANTIC_VERSION}.tar.gz"
E_ANTIC_SHA256=7546119231ba8fb6b3c9daf00c7979c1ec6f7692a75a8ee2ae6ac1abc5bddab4

CONFIGURE_FLAGS="--prefix=${PREFIX} --without-benchmark"

echo "Installing E-ANTIC..."

mkdir -p ${NMZ_OPT_DIR}/E-ANTIC_source/
cd ${NMZ_OPT_DIR}/E-ANTIC_source

../../download.sh ${E_ANTIC_URL} ${E_ANTIC_SHA256}
if [ ! -d e-antic-${E_ANTIC_VERSION} ]; then
    tar -xvf e-antic-${E_ANTIC_VERSION}.tar.gz
fi

cd e-antic-${E_ANTIC_VERSION}/libeantic

if [ ! -f config.status ]; then
    ./configure ${CONFIGURE_FLAGS}  CFLAGS="${CFLAGS} -I${PREFIX}/include" \
              CPPFLAGS="${CPPFLAGS}" LDFLAGS="${LDFLAGS}"
fi
make -j4
make install
