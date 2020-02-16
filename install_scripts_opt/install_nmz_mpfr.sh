#!/usr/bin/env bash

set -e

source $(dirname "$0")/common.sh

echo "Wauwau"
echo ${PREFIX}

CONFIGURE_FLAGS="--prefix=${PREFIX}"
if [ "$GMP_INSTALLDIR" != "" ]; then
    CONFIGURE_FLAGS="${CONFIGURE_FLAGS} --with-gmp=$GMP_INSTALLDIR"
fi

echo "MPFR flags"
echo $CONFIGURE_FLAGS

## script for the installation of MPFR (needed for Flint)

MPFR_VERSION="4.0.2"
MPFR_URL="https://ftp.gnu.org/gnu/mpfr/mpfr-${MPFR_VERSION}.tar.gz"
MPFR_SHA256=ae26cace63a498f07047a784cd3b0e4d010b44d2b193bab82af693de57a19a78

echo "Installing MPFR..."

# download & extract
mkdir -p ${NMZ_OPT_DIR}/MPFR_source/
cd ${NMZ_OPT_DIR}/MPFR_source
../../download.sh ${MPFR_URL} ${MPFR_SHA256}
if [ ! -d mpfr-${MPFR_VERSION} ]; then
    tar -xvf mpfr-${MPFR_VERSION}.tar.gz
fi

# configure & compile
cd mpfr-${MPFR_VERSION}
if [ ! -f config.status ]; then
    echo "Vor mpdfr configure"
    echo $CONFIGURE_FLAGS
    ./configure ${CONFIGURE_FLAGS}
fi

cat config.log

echo "Vor mpdfr make"

make -j4

echo "Vor mpdfr install"

make install

echo "Vor mpdfr done"
