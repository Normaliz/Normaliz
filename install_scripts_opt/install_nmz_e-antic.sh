#!/usr/bin/env bash

set -e

source $(dirname "$0")/common.sh

CONFIGURE_FLAGS="--prefix=${PREFIX}"
if [ "x$NO_OPENMP" != x ]; then
    CONFIGURE_FLAGS="${CONFIGURE_FLAGS} --disable-openmp"
fi

if [ "$GMP_INSTALLDIR" != "" ]; then
    CPPFLAGS="${CPFFLAGS} -I${GMP_INSTALLDIR}/include"
    LDFLAGS="${LDFLAGS} -L${GMP_INSTALLDIR}/lib"
fi

## script for the installation of e-antic for the use in libnormaliz
E_ANTIC_VERSION=0.1.3b0
E_ANTIC_URL="http://www.labri.fr/perso/vdelecro/e-antic/e-antic-${E_ANTIC_VERSION}.tar.gz"
E_ANTIC_SHA256=7bfe4aa926303b87a58a962535793bfadd9bbf0f7617e7be82e0f77e3351438e

echo "Installing E-ANTIC..."

# download & extract
mkdir -p ${NMZ_OPT_DIR}/E-ANTIC_source/
cd ${NMZ_OPT_DIR}/E-ANTIC_source

../../download.sh ${E_ANTIC_URL} ${E_ANTIC_SHA256}
if [ ! -d e-antic-${E_ANTIC_VERSION} ]; then
    tar -xvf e-antic-${E_ANTIC_VERSION}.tar.gz
fi

# configure & compile
cd e-antic-${E_ANTIC_VERSION}
if [ ! -f config.status ]; then
    ./configure ${CONFIGURE_FLAGS} CPPFLAGS="${CPPFLAGS} -fPIC"
# --enable-flint-devel ## for Flint development version
fi
make -j4
make install
