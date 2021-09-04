#!/usr/bin/env bash

set -e

echo "::group::nauty"

source $(dirname "$0")/common.sh

./install_scripts_opt/install_nmz_hash-library.sh

if [ "$GMP_INSTALLDIR" != "" ]; then
    CPPFLAGS="${CPPFLAGS} -I${GMP_INSTALLDIR}/include"
    LDFLAGS="${LDFLAGS} -L${GMP_INSTALLDIR}/lib"
fi

##  script for the installation of nauty
## as far as needed by libnormaliz

NAUTY_VERSION="27r2"
# NAUTY_URL="https://pallini.di.uniroma1.it/nauty${NAUTY_VERSION}.tar.gz"
NAUTY_URL="https://users.cecs.anu.edu.au/~bdm/nauty/nauty${NAUTY_VERSION}.tar.gz"
NAUTY_SHA256=fc434729c833d9bb7053c8def10e72eaede03487d1afa50568ce2972b0337741

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
./configure --enable-tls
make all -j4 CFLAGS="-fPIC -O3"
mkdir -p ${PREFIX}/include/nauty
cp nauty.h ${PREFIX}/include/nauty
# mkdir -p ${PREFIX}/lib ## in common.sh
cp nauty.a ${PREFIX}/lib/libnauty.a

echo "nauty installed"
