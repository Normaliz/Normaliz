#!/usr/bin/env bash

set -e

source $(dirname "$0")/common.sh

echo "Prefix rb"
echo ${PREFIX}

echo "Start Arb"
ls ${PREFIX}/lib

CONFIGURE_FLAGS="--prefix=${PREFIX}"
if [ "$GMP_INSTALLDIR" != "" ]; then
    CONFIGURE_FLAGS="${CONFIGURE_FLAGS} --with-gmp=$GMP_INSTALLDIR"
fi

# unfortunately arb does not support standard env vars like CPPFLAGS, LDFLAGS;
# but luckily it has some vars of its own that we can use for the same purpose
export EXTRA_INC_DIRS="${PREFIX}/include"
export EXTRA_LIB_DIRS="${PREFIX}/lib"

## script for the installation of ARB for the use in libnormaliz

ARB_VERSION="2.16.0"
ARB_URL="https://github.com/fredrik-johansson/arb/archive/${ARB_VERSION}.tar.gz"
ARB_SHA256=77464be4d34a511bb004457f862fec857ff934b0ed58d56d6f52d76ebadd4daf

echo "Installing ARB..."

# download & extract
mkdir -p ${NMZ_OPT_DIR}/ARB_source/
cd ${NMZ_OPT_DIR}/ARB_source
../../download.sh ${ARB_URL} ${ARB_SHA256} arb-${ARB_VERSION}.tar.gz
if [ ! -d arb-${ARB_VERSION} ]; then
    tar -xvf arb-${ARB_VERSION}.tar.gz
fi

# configure & compile
cd arb-${ARB_VERSION}
if [ ! -f Makefile ]; then
    ./configure ${CONFIGURE_FLAGS}
fi
make -j4 verbose
make install

echo "End Arb"
ls ${PREFIX}/lib
