#!/usr/bin/env bash

set -e

if [ "$GMP_INSTALLDIR" != "" ]; then
    CPPFLAGS="${CPFFLAGS} -I${GMP_INSTALLDIR}/include"
    LDFLAGS="${LDFLAGS} -L${GMP_INSTALLDIR}/lib"
fi

source $(dirname "$0")/common.sh

##  script for the installation of CoCoALib
## as far as needed by libnormaliz

COCOA_VERSION="0.99601"
COCOA_URL="http://cocoa.dima.unige.it/cocoalib/tgz/CoCoALib-${COCOA_VERSION}.tgz"
COCOA_SHA256=caf37f71398b9715be262e434f04a218db05cfa58e08bce954626d7f4ffd6b75

echo "Installing CoCoA..."

mkdir -p ${NMZ_OPT_DIR}/CoCoA_source/
cd ${NMZ_OPT_DIR}/CoCoA_source/
../../download.sh ${COCOA_URL} ${COCOA_SHA256}
if [ ! -d CoCoALib-${COCOA_VERSION} ]; then
    tar xvf CoCoALib-${COCOA_VERSION}.tgz
fi
cd CoCoALib-${COCOA_VERSION}
if [ ! -f configuration/autoconf.mk ]; then
    ./configure --threadsafe-hack --no-boost
fi
make library -j4
mkdir -p ${PREFIX}/include/CoCoA
cp include/CoCoA/*.H ${PREFIX}/include/CoCoA
mkdir -p ${PREFIX}/lib
cp lib/libcocoa.a ${PREFIX}/lib
