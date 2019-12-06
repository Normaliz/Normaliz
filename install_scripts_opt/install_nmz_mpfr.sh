#!/usr/bin/env bash

set -e

WITH_GMP=""
if [ "$GMP_INSTALLDIR" != "" ]; then
  WITH_GMP="--with-gmp=$GMP_INSTALLDIR"
fi

source $(dirname "$0")/common.sh

## script for the installation of MPFR (needed for Flint)

MPFR_VERSION="4.0.2"
MPFR_URL="https://ftp.gnu.org/gnu/mpfr/mpfr-${MPFR_VERSION}.tar.gz"
MPFR_SHA256=ae26cace63a498f07047a784cd3b0e4d010b44d2b193bab82af693de57a19a78

echo "Installing MPFR..."

mkdir -p ${NMZ_OPT_DIR}/MPFR_source/
cd ${NMZ_OPT_DIR}/MPFR_source
../../download.sh ${MPFR_URL} ${MPFR_SHA256}
if [ ! -d mpfr-${MPFR_VERSION} ]; then
    tar -xvf mpfr-${MPFR_VERSION}.tar.gz
fi
cd mpfr-${MPFR_VERSION}
if [ ! -f config.status ]; then
    ./configure --prefix=${PREFIX} $WITH_GMP
fi
make -j4
make install
