#!/usr/bin/env bash

set -e

WITH_GMP=""
if [ "$GMP_INSTALLDIR" != "" ]; then
  WITH_GMP="--with-gmp=$GMP_INSTALLDIR"
fi

if [ "x$NMZ_OPT_DIR" = x ]; then 
    export NMZ_OPT_DIR=${PWD}/nmz_opt_lib
    mkdir -p ${NMZ_OPT_DIR}
fi

if [ "x$NMZ_COMPILER" != x ]; then
    export CXX=$NMZ_COMPILER
elif [[ $OSTYPE == darwin* ]]; then
    export CXX=clang++
    export PATH="`brew --prefix`/opt/llvm/bin/:$PATH"
    export LDFLAGS="${LDFLAGS} -L`brew --prefix`/opt/llvm/lib"
fi

## first install mpfr
$(dirname "$0")/install_nmz_mpfr.sh

## script for the installation of Flint for the use in libnormaliz

FLINT_VERSION="2.5.2"
FLINT_URL="http://www.flintlib.org/flint-${FLINT_VERSION}.tar.gz"
FLINT_SHA256=cbf1fe0034533c53c5c41761017065f85207a1b770483e98b2392315f6575e87

if [ "x$NMZ_PREFIX" != x ]; then
    mkdir -p ${NMZ_PREFIX}
    PREFIX=${NMZ_PREFIX}
else
    PREFIX=${PWD}/local
fi

echo "Installing FLINT..."

mkdir -p ${NMZ_OPT_DIR}/Flint_source/
cd ${NMZ_OPT_DIR}/Flint_source
../../download.sh ${FLINT_URL} ${FLINT_SHA256}
if [ ! -d flint-${FLINT_VERSION} ]; then
    tar -xvf flint-${FLINT_VERSION}.tar.gz
fi
cd flint-${FLINT_VERSION}
if [ ! -f Makefile ]; then
    ./configure --prefix=${PREFIX} --with-mpfr=${PREFIX} $WITH_GMP $EXTRA_FLINT_FLAGS
fi
# patch to avoid PIE clash in Ubuntu >= 16-10
if [[ $OSTYPE == "linux-gnu" ]]; then
    sed -i s/"-Wl,"// Makefile.subdirs
fi
make -j4 verbose
make install
