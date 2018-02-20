#!/usr/bin/env bash -e

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

##  script for the installation of CoCoALib
## as far as needed by libnormaliz

COCOA_VERSION="0.99550"

INSTALLDIR=${NMZ_OPT_DIR}

mkdir -p ${NMZ_OPT_DIR}/CoCoA_source/
cd ${NMZ_OPT_DIR}/CoCoA_source/
curl -O http://cocoa.dima.unige.it/cocoalib/tgz/CoCoALib-${COCOA_VERSION}.tgz
tar xvf CoCoALib-${COCOA_VERSION}.tgz
cd CoCoALib-${COCOA_VERSION}
./configure --threadsafe-hack --no-boost $WITH_GMP
make library -j4
mkdir -p  ${INSTALLDIR}/include
mkdir  -p ${INSTALLDIR}/include/CoCoA
cp include/CoCoA/*.H  ${INSTALLDIR}/include/CoCoA
mkdir  -p ${INSTALLDIR}/lib
cp lib/libcocoa.a ${INSTALLDIR}/lib
