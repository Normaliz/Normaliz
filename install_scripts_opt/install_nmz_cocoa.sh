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

##  script for the installation of CoCoALib
## as far as needed by libnormaliz

COCOA_VERSION="0.99601"
COCOA_URL="http://cocoa.dima.unige.it/cocoalib/tgz/CoCoALib-${COCOA_VERSION}.tgz"
COCOA_SHA256=caf37f71398b9715be262e434f04a218db05cfa58e08bce954626d7f4ffd6b75

#INSTALLDIR=${NMZ_OPT_DIR}
if [ "x$NMZ_PREFIX" != x ]; then
    mkdir -p ${NMZ_PREFIX}
    INSTALLDIR=${NMZ_PREFIX}
else
    INSTALLDIR=${PWD}/local
fi


echo "Installing CoCoA..."

mkdir -p ${NMZ_OPT_DIR}/CoCoA_source/
cd ${NMZ_OPT_DIR}/CoCoA_source/
../../download.sh ${COCOA_URL} ${COCOA_SHA256}
if [ ! -d CoCoALib-${COCOA_VERSION} ]; then
    tar xvf CoCoALib-${COCOA_VERSION}.tgz
fi
cd CoCoALib-${COCOA_VERSION}
if [ ! -f configuration/autoconf.mk ]; then
    ./configure --threadsafe-hack --no-boost
fi
make library -j4
mkdir -p  ${INSTALLDIR}/include
mkdir  -p ${INSTALLDIR}/include/CoCoA
cp include/CoCoA/*.H  ${INSTALLDIR}/include/CoCoA
mkdir  -p ${INSTALLDIR}/lib
cp lib/libcocoa.a ${INSTALLDIR}/lib
