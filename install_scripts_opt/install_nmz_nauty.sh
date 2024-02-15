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

NAUTY_VERSION="2_8_8"
# NAUTY_URL="https://pallini.di.uniroma1.it/nauty${NAUTY_VERSION}.tar.gz"
NAUTY_URL="https://users.cecs.anu.edu.au/~bdm/nauty/nauty${NAUTY_VERSION}.tar.gz"
NAUTY_SHA256=159d2156810a6bb240410cd61eb641add85088d9f15c888cdaa37b8681f929ce

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
if [ "$OSTYPE" != "msys" ]; then
	cp nauty.a ${PREFIX}/lib/libnauty.a ## WORDSIZE = 64
else
	cp nautyW.a ${PREFIX}/lib/libnauty.a ## WORDSIZE = 32
fi


echo "nauty installed"
