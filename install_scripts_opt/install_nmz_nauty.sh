#!/usr/bin/env bash

set -e

source $(dirname "$0")/common.sh

if [ "$GMP_INSTALLDIR" != "" ]; then
    CPPFLAGS="${CPFFLAGS} -I${GMP_INSTALLDIR}/include"
    LDFLAGS="${LDFLAGS} -L${GMP_INSTALLDIR}/lib"
fi

##  script for the installation of nauty
## as far as needed by libnormaliz

NAUTY_VERSION="27rc3"
NAUTY_URL="http://pallini.di.uniroma1.it/nauty${NAUTY_VERSION}.tar.gz"
NAUTY_SHA256=3293929ce05efc338dc9c21dc3db41751c55d3ae9a78dfb31fa3195839f4e0c7

echo "Installing nauty..."

# download & extract
mkdir -p ${NMZ_OPT_DIR}/Nauty_source/
cd ${NMZ_OPT_DIR}/Nauty_source/
../../download.sh ${NAUTY_URL} ${NAUTY_SHA256}
if [ ! -d nauty${NAUTY_VERSION} ]; then
    tar xvf nauty${NAUTY_VERSION}.tar.gz
fi
cd nauty${NAUTY_VERSION}

# configure & compile
./configure
make all -j4 CFLAGS="-fPIC -O3"
mkdir -p ${PREFIX}/include/nauty
cp nauty.h ${PREFIX}/include/nauty
mkdir -p ${PREFIX}/lib
cp nauty.a ${PREFIX}/lib/libnauty.a
