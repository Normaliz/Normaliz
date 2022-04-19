#!/usr/bin/env bash

## script for the installation of ARB for the use in libnormaliz

set -e

echo "::group::arb"

source $(dirname "$0")/common.sh

CONFIGURE_FLAGS="--prefix=${PREFIX}"
if [ "$OSTYPE" != "msys" ]; then
	CONFIGURE_FLAGS="${CONFIGURE_FLAGS} --with-mpfr=${PREFIX} --with-flint=${PREFIX}"
else
	CONFIGURE_FLAGS="${CONFIGURE_FLAGS} --with-gmp=${MSYS_STANDARD_LOC} --with-mpfr=${MSYS_STANDARD_LOC} --with-flint=${MSYS_STANDARD_LOC}"
fi

if [ "$GMP_INSTALLDIR" != "" ]; then
    CONFIGURE_FLAGS="${CONFIGURE_FLAGS} --with-gmp=${GMP_INSTALLDIR}"
fi

ARB_VERSION="2.22.0"
ARB_URL="https://github.com/fredrik-johansson/arb/archive/${ARB_VERSION}.tar.gz"
ARB_SHA256=3e40ab8cf61c0cd63d5901064d73eaa2d04727bbdc6eebb1727997958a14f24d

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
if [ "$OSTYPE" != "msys" ]; then
	./configure ${CONFIGURE_FLAGS}
	make -j4
	make install
else
	# must make shared and static in separate builds
	./configure ${CONFIGURE_FLAGS} --disable-static
	make -j4
	make install
	./configure ${CONFIGURE_FLAGS} --disable-shared
	make -j4
	make install
fi
