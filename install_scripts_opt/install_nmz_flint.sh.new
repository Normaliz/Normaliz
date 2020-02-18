#!/usr/bin/env bash

set -e

source $(dirname "$0")/common.sh

CONFIGURE_FLAGS="--prefix=${PREFIX}"
if [ "$GMP_INSTALLDIR" != "" ]; then
    CONFIGURE_FLAGS="${CONFIGURE_FLAGS} --with-gmp=$GMP_INSTALLDIR"
fi

# unfortunately flint does not support standard env vars like CPPFLAGS, LDFLAGS;
# but luckily it has some vars of its own that we can use for the same purpose
export EXTRA_INC_DIRS="${PREFIX}/include"
export EXTRA_LIB_DIRS="${PREFIX}/lib"

## script for the installation of Flint for the use in libnormaliz

FLINT_VERSION="2.5.2"
FLINT_URL="http://www.flintlib.org/flint-${FLINT_VERSION}.tar.gz"
FLINT_SHA256=cbf1fe0034533c53c5c41761017065f85207a1b770483e98b2392315f6575e87

echo "Installing Flint..."

# download & extract
mkdir -p ${NMZ_OPT_DIR}/Flint_source/
cd ${NMZ_OPT_DIR}/Flint_source
../../download.sh ${FLINT_URL} ${FLINT_SHA256}
if [ ! -d flint-${FLINT_VERSION} ]; then
    tar -xvf flint-${FLINT_VERSION}.tar.gz
fi

# configure & compile
cd flint-${FLINT_VERSION}
if [ ! -f Makefile ]; then
    ./configure ${CONFIGURE_FLAGS} $EXTRA_FLINT_FLAGS
fi
# patch to avoid PIE clash in Ubuntu >= 16-10
if [[ $OSTYPE == "linux-gnu" ]]; then
    sed -i s/"-Wl,"// Makefile.subdirs
fi
make -j4
make install

echo "Flint installed"

