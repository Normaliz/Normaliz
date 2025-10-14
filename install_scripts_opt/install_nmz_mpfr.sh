#!/usr/bin/env bash

set -e

echo "::group::mpfr"

source $(dirname "$0")/common.sh

CONFIGURE_FLAGS="--prefix=${PREFIX}"

if [ "$OSTYPE" != "msys" ]; then
	CONFIGURE_FLAGS="${CONFIGURE_FLAGS} --with-mpfr=${PREFIX}"
else # only static here, we take shared from MSYS repository
	CONFIGURE_FLAGS="${CONFIGURE_FLAGS} --with-gmp=${MSYS_STANDARD_LOC} --disable-shared"
fi

if [ "$GMP_INSTALLDIR" != "" ]; then
    CONFIGURE_FLAGS="${CONFIGURE_FLAGS} --with-gmp=$GMP_INSTALLDIR"
fi


echo "MPFR flags"
echo $CONFIGURE_FLAGS

## script for the installation of MPFR (needed for Flint)

MPFR_VERSION="4.2.2"
MPFR_URL="https://ftp.gnu.org/gnu/mpfr/mpfr-${MPFR_VERSION}.tar.xz"
MPFR_SHA256=b67ba0383ef7e8a8563734e2e889ef5ec3c3b898a01d00fa0a6869ad81c6ce01

echo "Installing MPFR..."

# download & extract
mkdir -p ${NMZ_OPT_DIR}/MPFR_source/
cd ${NMZ_OPT_DIR}/MPFR_source
../../download.sh ${MPFR_URL} ${MPFR_SHA256}
if [ ! -d mpfr-${MPFR_VERSION} ]; then
    tar -xvf mpfr-${MPFR_VERSION}.tar.xz
fi

# configure & compile
cd mpfr-${MPFR_VERSION}
if [ ! -f config.status ]; then
    echo "Vor mpfr configure"
    echo $CONFIGURE_FLAGS
    ./configure ${CONFIGURE_FLAGS}
fi

make -j4
make install

echo "MPFR installed"
