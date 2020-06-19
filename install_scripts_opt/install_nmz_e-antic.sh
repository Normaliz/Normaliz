#!/usr/bin/env bash

set -e

if [ "$GMP_INSTALLDIR" != "" ]; then
    CPPFLAGS="${CPFFLAGS} -I${GMP_INSTALLDIR}/include"
    LDFLAGS="${LDFLAGS} -L${GMP_INSTALLDIR}/lib"
fi

source $(dirname "$0")/common.sh

## script for the installation of e-antic for the use in libnormaliz
E_ANTIC_VERSION=0.1.7
E_ANTIC_URL="http://www.labri.fr/perso/vdelecro/e-antic/e-antic-${E_ANTIC_VERSION}.tar.gz"
E_ANTIC_SHA256=031a7702296b41a76b17ea9fbca184ac6ffae539d7ed0bdd987c0db1ab825f0a

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
