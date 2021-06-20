#!/usr/bin/env bash

set -e

echo "::group::antic"

source $(dirname "$0")/common.sh

CONFIGURE_FLAGS="--prefix=${PREFIX} --with-mpfr=${PREFIX} --with-flint=${PREFIX}"
if [ "$GMP_INSTALLDIR" != "" ]; then
    CONFIGURE_FLAGS="${CONFIGURE_FLAGS} --with-gmp=${GMP_INSTALLDIR}"
else
    CONFIGURE_FLAGS="${CONFIGURE_FLAGS} --with-gmp=${PREFIX}"
fi

## script for the installation of e-antic for the use in libnormaliz
ANTIC_VERSION=0.2.4
ANTIC_URL="https://github.com/wbhart/antic/archive/v${ANTIC_VERSION}.tar.gz"
ANTIC_SHA256=517d53633ff9c6348549dc6968567051b2161098d2bc395cb40ecc41e24312c6

mkdir -p ${NMZ_OPT_DIR}/ANTIC_source/
cd ${NMZ_OPT_DIR}/ANTIC_source

../../download.sh ${ANTIC_URL} ${ANTIC_SHA256}
if [ ! -d antic-${ANTIC_VERSION} ]; then
    tar -xvf v${ANTIC_VERSION}.tar.gz
fi

cd antic-${ANTIC_VERSION}/

./configure ${CONFIGURE_FLAGS}

make -j4
make install

ls -l ${PREFIX}/lib
