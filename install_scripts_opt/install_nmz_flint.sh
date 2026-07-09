#!/usr/bin/env bash

set -e

echo "::group::flint"

source $(dirname "$0")/common.sh

CONFIGURE_FLAGS="--prefix=${PREFIX}"

if [ "$OSTYPE" != "msys" ]; then
	CONFIGURE_FLAGS="${CONFIGURE_FLAGS} --with-mpfr=${PREFIX} --enable-static"
else # only static here, we take shared from MSYS repository
	CONFIGURE_FLAGS="${CONFIGURE_FLAGS} --with-gmp=${MSYS_STANDARD_LOC} --with-mpfr=${MSYS_STANDARD_LOC} --disable-shared"
fi

if [ "$GMP_INSTALLDIR" != "" ]; then
    CONFIGURE_FLAGS="${CONFIGURE_FLAGS} --with-gmp=${GMP_INSTALLDIR}"
fi

FLINT_CFLAGS="${FLINT_CFLAGS:--O2 -g}"

## script for the installation of Flint for the use in libnormaliz

FLINT_VERSION="3.3.1"
FLINT_URL="https://flintlib.org/download/flint-${FLINT_VERSION}.tar.gz"
FLINT_SHA256=64d70e513076cfa971e0410b58c1da5d35112913e9a56b44e2c681b459d3eafb

echo "Installing FLINT..."

mkdir -p ${NMZ_OPT_DIR}/Flint_source/
cd ${NMZ_OPT_DIR}/Flint_source
../../download.sh ${FLINT_URL} ${FLINT_SHA256}
if [ ! -d flint-${FLINT_VERSION} ]; then
    tar -xvf flint-${FLINT_VERSION}.tar.gz
fi
cd flint-${FLINT_VERSION}
if [ ! -f Makefile ]; then
    ./configure ${CONFIGURE_FLAGS} CFLAGS="${FLINT_CFLAGS}"
fi
if grep -E '(^|[^[:alnum:]_])-m(arch|tune|avx)' Makefile; then
    echo "error: FLINT configure emitted CPU-specific compiler flags"
    exit 1
fi
# patch to avoid PIE clash in Ubuntu >= 16-10
## if [[ $OSTYPE == "linux-gnu" ]]; then
## sed -i s/"-Wl,"// Makefile.subdirs
## fi
# make -j4 # verbose
make install -j8
