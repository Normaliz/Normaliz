#!/usr/bin/env bash

set -e

WITH_GMP=""
if [ "$GMP_INSTALLDIR" != "" ]; then
  WITH_GMP="--with-libgmp=$GMP_INSTALLDIR/lib/libgmp.a"
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

#INSTALLDIR=${NMZ_OPT_DIR}
if [ "x$NMZ_PREFIX" != x ]; then
    mkdir -p ${NMZ_PREFIX}
    INSTALLDIR=${NMZ_PREFIX}
else
    INSTALLDIR=${PWD}/local
fi

##  script for the installation of nauty
## as far as needed by libnormaliz

NAUTY_VERSION="27rc2"


echo "Installing nauty..."

mkdir -p ${NMZ_OPT_DIR}/Nauty_source/
cd ${NMZ_OPT_DIR}/Nauty_source/
if [ ! -d nauty${NAUTY_VERSION} ]; then
    wget http://pallini.di.uniroma1.it/nauty${NAUTY_VERSION}.tar.gz
    tar xvf nauty${NAUTY_VERSION}.tar.gz
fi
cd nauty${NAUTY_VERSION}

    ./configure

make all -j4 CFLAGS="-fPIC -O3"
mkdir -p  ${INSTALLDIR}/include
mkdir  -p ${INSTALLDIR}/include/nauty
cp nauty.h  ${INSTALLDIR}/include/nauty
mkdir  -p ${INSTALLDIR}/lib
cp nauty.a ${INSTALLDIR}/lib/libnauty.a
