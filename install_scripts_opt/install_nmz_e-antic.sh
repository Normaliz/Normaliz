#!/usr/bin/env bash

set -e

if [ "$GMP_INSTALLDIR" != "" ]; then
    CPPFLAGS="${CPFFLAGS} -I${GMP_INSTALLDIR}/include"
    LDFLAGS="${LDFLAGS} -L${GMP_INSTALLDIR}/lib"
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

## script for the installation of e-antic for the use in libnormaliz
E_ANTIC_VERSION=0.1.3b0
E_ANTIC_URL="http://www.labri.fr/perso/vdelecro/e-antic/e-antic-${E_ANTIC_VERSION}.tar.gz"
E_ANTIC_SHA256=7bfe4aa926303b87a58a962535793bfadd9bbf0f7617e7be82e0f77e3351438e

if [ "x$NMZ_PREFIX" != x ]; then
    mkdir -p ${NMZ_PREFIX}
    PREFIX=${NMZ_PREFIX}
else
    PREFIX=${PWD}/local
fi

if [ "x$NO_OPENMP" != x ]; then
    export BLOCK_OPENMP="--disable-openmp"
fi

echo "Installing E-ANTIC..."

mkdir -p ${NMZ_OPT_DIR}/E-ANTIC_source/
cd ${NMZ_OPT_DIR}/E-ANTIC_source

../../download.sh ${E_ANTIC_URL} ${E_ANTIC_SHA256}
if [ ! -d e-antic-${E_ANTIC_VERSION} ]; then
    tar -xvf e-antic-${E_ANTIC_VERSION}.tar.gz

fi

cd e-antic-${E_ANTIC_VERSION}

if [ ! -f config.status ]; then
    ./configure --prefix=${PREFIX}  ${BLOCK_OPENMP} \
              CFLAGS="${CFLAGS} -I${PREFIX}/include" \
              CPPFLAGS="${CPPFLAGS} -I${PREFIX}/include -fPIC" \
              LDFLAGS="${LDFLAGS} -L/${PREFIX}/lib"
# --enable-flint-devel ## for Flint development version
fi
make -j4
make install
