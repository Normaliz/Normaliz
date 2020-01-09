#!/usr/bin/env bash

set -e

source $(dirname "$0")/common.sh

if [ "$GMP_INSTALLDIR" != "" ]; then
    CPPFLAGS="${CPFFLAGS} -I${GMP_INSTALLDIR}/include"
    LDFLAGS="${LDFLAGS} -L${GMP_INSTALLDIR}/lib"
fi

##  script for the installation of nauty
## as far as needed by libnormaliz

NAUTY_VERSION="27rc2"
NAUTY_URL="http://pallini.di.uniroma1.it/nauty${NAUTY_VERSION}.tar.gz"
NAUTY_SHA256=fe5b893579d84736a6e1383e0dfe95f9df55bff6ad3eaafd44660e08745ab32f

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
