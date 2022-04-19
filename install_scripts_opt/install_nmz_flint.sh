#!/usr/bin/env bash

set -e

echo "::group::flint"

source $(dirname "$0")/common.sh

CONFIGURE_FLAGS="--prefix=${PREFIX}"

if [ "$OSTYPE" != "msys" ]; then
	CONFIGURE_FLAGS="${CONFIGURE_FLAGS} --with-mpfr=${PREFIX}"
else # only static here, we take shared from MSYS repository
	CONFIGURE_FLAGS="${CONFIGURE_FLAGS} --with-gmp=${MSYS_STANDARD_LOC} --with-mpfr=${MSYS_STANDARD_LOC} --disable-shared"
fi

if [ "$GMP_INSTALLDIR" != "" ]; then
    CONFIGURE_FLAGS="${CONFIGURE_FLAGS} --with-gmp=${GMP_INSTALLDIR}"
fi

## script for the installation of Flint for the use in libnormaliz

FLINT_VERSION="2.8.4"
FLINT_URL="https://flintlib.org/flint-${FLINT_VERSION}.tar.gz"
FLINT_SHA256=61df92ea8c8e9dc692d46c71d7f50aaa09a33d4ba08d02a1784730a445e5e4be

echo "Installing FLINT..."

mkdir -p ${NMZ_OPT_DIR}/Flint_source/
cd ${NMZ_OPT_DIR}/Flint_source
../../download.sh ${FLINT_URL} ${FLINT_SHA256}
if [ ! -d flint-${FLINT_VERSION} ]; then
    tar -xvf flint-${FLINT_VERSION}.tar.gz
fi
cd flint-${FLINT_VERSION}
if [ ! -f Makefile ]; then
    ./configure ${CONFIGURE_FLAGS}
fi
# patch to avoid PIE clash in Ubuntu >= 16-10
if [[ $OSTYPE == "linux-gnu" ]]; then
    sed -i s/"-Wl,"// Makefile.subdirs
fi
# make -j4 # verbose
make install -j8
