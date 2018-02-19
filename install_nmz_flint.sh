#!/bin/sh -e

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
    export LDFLAGS="-L`brew --prefix`/opt/llvm/lib"
fi

## script for the installation of Flint for the use in libnormaliz
## including the installation of MPFR (needed for Flint)

FLINT_VERSION="2.5.2"
MPFR_VERSION="4.0.0"

PREFIX=${NMZ_OPT_DIR}

mkdir -p ${NMZ_OPT_DIR}/MPFR_source/
cd ${NMZ_OPT_DIR}/MPFR_source
curl -O http://www.mpfr.org/mpfr-current/mpfr-${MPFR_VERSION}.tar.gz 
tar -xvf mpfr-${MPFR_VERSION}.tar.gz
cd mpfr-${MPFR_VERSION}
./configure --prefix=${PREFIX} $WITH_GMP
make -j4
make install

mkdir -p ${NMZ_OPT_DIR}/Flint_source/
cd ${NMZ_OPT_DIR}/Flint_source
curl -O http://www.flintlib.org/flint-${FLINT_VERSION}.tar.gz
tar -xvf flint-${FLINT_VERSION}.tar.gz
cd flint-${FLINT_VERSION}
./configure --prefix=${PREFIX} --with-mpfr=${PREFIX} $WITH_GMP $EXTRA_FLINT_FLAGS
make -j4
make install
