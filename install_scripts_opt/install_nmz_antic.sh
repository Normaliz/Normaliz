#!/usr/bin/env bash

set -e

source $(dirname "$0")/common.sh

if [ "$GMP_INSTALLDIR" != "" ]; then
    CONFIGURE="--with-gmp=${GMP_INSTALLDIR}"
else
    CONFIGURE="--with-gmp=${PREFIX}"
fi

## script for the installation of e-antic for the use in libnormaliz
ANTIC_VERSION=0.2.1
ANTIC_URL="https://github.com/wbhart/antic/archive/v${ANTIC_VERSION}.tar.gz"
ANTIC_SHA256=b44dfcaa93db9e9a9fb3e47ce861e9e74de6d53e53754dd821b80b3e02e27607

CONFIGURE_FLAGS="--prefix=${PREFIX} ${CONFIGURE} --with-mpfr=${PREFIX} --with-flint=${PREFIX}"

echo "Installing ANTIC..."

mkdir -p ${NMZ_OPT_DIR}/ANTIC_source/
cd ${NMZ_OPT_DIR}/ANTIC_source

../../download.sh ${ANTIC_URL} ${ANTIC_SHA256}
if [ ! -d antic-${ANTIC_VERSION} ]; then
    tar -xvf antic-${ANTIC_VERSION}.tar.gz
fi

cd antic-${ANTIC_VERSION}/

./configure ${CONFIGURE_FLAGS} CFLAGS="${CFLAGS} -I${PREFIX}/include" CPPFLAGS="${CPPFLAGS}" LDFLAGS="${LDFLAGS}"

make -j4
make install

