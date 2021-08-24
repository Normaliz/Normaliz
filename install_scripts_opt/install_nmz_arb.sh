#!/usr/bin/env bash

set -e

echo "::group::arb"

source $(dirname "$0")/common.sh

CONFIGURE_FLAGS="--prefix=${PREFIX} --with-mpfr=${PREFIX} --with-flint=${PREFIX} ${EXTRA_ARB_FLAGS}"
if [ "$GMP_INSTALLDIR" != "" ]; then
    CONFIGURE_FLAGS="${CONFIGURE_FLAGS} --with-gmp=${GMP_INSTALLDIR}"
fi

## script for the installation of ARB for the use in libnormaliz

ARB_VERSION="2.20.0"
ARB_URL="https://github.com/fredrik-johansson/arb/archive/${ARB_VERSION}.tar.gz"
ARB_SHA256=d2f186b10590c622c11d1ca190c01c3da08bac9bc04e84cb591534b917faffe7

echo "Installing ARB..."

mkdir -p ${NMZ_OPT_DIR}/ARB_source/
cd ${NMZ_OPT_DIR}/ARB_source
# ../../download.sh ${ARB_URL} ${ARB_SHA256} arb-${ARB_VERSION}.tar.gz
../../download.sh ${ARB_URL} ${ARB_SHA256} arb-${ARB_VERSION}.tar.gz
if [ ! -d arb-${ARB_VERSION} ]; then
    tar -xvf arb-${ARB_VERSION}.tar.gz
fi
cd arb-${ARB_VERSION}
# (In particular on Mac OS X, make sure that our version of MPFR comes
# first in the -L search path, not the one from LLVM or elsewhere.
# ARB's configure puts it last.)
## export LDFLAGS="-L${NMZ_OPT_DIR}/lib ${LDFLAGS}"
# export LDFLAGS="-L${PREFIX}/lib ${LDFLAGS}"
if [ ! -f Makefile ]; then
    ./configure ${CONFIGURE_FLAGS}
fi
make -j4 # verbose
make install
