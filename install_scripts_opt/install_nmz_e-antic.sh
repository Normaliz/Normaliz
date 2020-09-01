#!/usr/bin/env bash

set -e

if [ "$GMP_INSTALLDIR" != "" ]; then
    CPPFLAGS="${CPFFLAGS} -I${GMP_INSTALLDIR}/include"
    LDFLAGS="${LDFLAGS} -L${GMP_INSTALLDIR}/lib"
fi

source $(dirname "$0")/common.sh

## script for the installation of e-antic for the use in libnormaliz
E_ANTIC_VERSION=0.1.8
E_ANTIC_URL="http://www.labri.fr/perso/vdelecro/e-antic/e-antic-${E_ANTIC_VERSION}.tar.gz"
E_ANTIC_SHA256=d8b6c18107756db86c7c4ca1593364e0ba9b4ad125910b45afd7ad5823d7100c

CONFIGURE_FLAGS="--prefix=${PREFIX}"

if [ "x$NO_OPENMP" != x ]; then
    CONFIGURE_FLAGS="${CONFIGURE_FLAGS} --disable-openmp"
fi

echo "Installing E-ANTIC..."

mkdir -p ${NMZ_OPT_DIR}/E-ANTIC_source/
cd ${NMZ_OPT_DIR}/E-ANTIC_source

../../download.sh ${E_ANTIC_URL} ${E_ANTIC_SHA256}
if [ ! -d e-antic-${E_ANTIC_VERSION} ]; then
    tar -xvf e-antic-${E_ANTIC_VERSION}.tar.gz

fi

cd e-antic-${E_ANTIC_VERSION}

# copy patch for Flint 2.6.1
cp ../../../install_scripts_opt/e-antic_patch/*.h e-antic
cp ../../../install_scripts_opt/e-antic_patch/*.c nf_elem

if [ ! -f config.status ]; then
    ./configure ${CONFIGURE_FLAGS}  CFLAGS="${CFLAGS} -I${PREFIX}/include" \
              CPPFLAGS="${CPPFLAGS}" LDFLAGS="${LDFLAGS}"
# --enable-flint-devel ## for Flint development version
fi
make -j4
make install
