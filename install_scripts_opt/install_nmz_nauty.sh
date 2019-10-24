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
NAUTY_URL="http://pallini.di.uniroma1.it/nauty${NAUTY_VERSION}.tar.gz"
NAUTY_SHA256=fe5b893579d84736a6e1383e0dfe95f9df55bff6ad3eaafd44660e08745ab32f


echo "Installing nauty..."

mkdir -p ${NMZ_OPT_DIR}/Nauty_source/
cd ${NMZ_OPT_DIR}/Nauty_source/
../../download.sh ${NAUTY_URL} ${NAUTY_SHA256}
if [ ! -d nauty${NAUTY_VERSION} ]; then
    tar xvf nauty${NAUTY_VERSION}.tar.gz
fi
cd nauty${NAUTY_VERSION}

./configure

make all -j4 CFLAGS="-fPIC -O3 -mpopcnt -march=native"
mkdir -p ${INSTALLDIR}/include
mkdir -p ${INSTALLDIR}/include/nauty
cp nauty.h ${INSTALLDIR}/include/nauty
mkdir -p ${INSTALLDIR}/lib
cp nauty.a ${INSTALLDIR}/lib/libnauty.a
