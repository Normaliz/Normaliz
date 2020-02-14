#!/usr/bin/env bash

set -e

source $(dirname "$0")/install_scripts_opt/common.sh

WITH_GMP=""
if [ "$GMP_INSTALLDIR" != "" ]; then
  WITH_GMP="--with-gmp=$GMP_INSTALLDIR"
fi

if [ ! -e configure ]; then
    ./bootstrap.sh
fi

if [ "x$NO_OPENMP" != x ]; then
    export BLOCK_OPENMP="--disable-openmp"
fi

echo "installing shared"

mkdir -p build
cd build

if [ "x$NMZ_EXTENDED_TESTS" != x ]; then
../configure --prefix="${PREFIX}"  CPPFLAGS="-I ${PREFIX}/include -DNMZ_EXTENDED_TESTS" LDFLAGS="-L${PREFIX}/lib/" $EXTRA_FLAGS $WITH_GMP ${BLOCK_OPENMP} --srcdir=..
else
../configure --prefix="${PREFIX}"  CPPFLAGS="-I ${PREFIX}/include" LDFLAGS="-L${PREFIX}/lib/" $EXTRA_FLAGS $WITH_GMP ${BLOCK_OPENMP} --srcdir=..
fi

make clean
make -j4
make install

# make distclean

cd ..
# rm -r build

cp -f ${PREFIX}/bin/* .
cp ${PREFIX}/lib/libnormaliz.a source/libnormaliz ## for compatibility with Makefile.classic

echo "Normaliz installation complete"
