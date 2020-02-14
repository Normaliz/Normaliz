#!/usr/bin/env bash

set -e

source $(dirname "$0")/install_scripts_opt/common.sh

if [ ! -e configure ]; then
    ./bootstrap.sh
fi

CONFIGURE_FLAGS="--prefix=${PREFIX}"
if [ "$GMP_INSTALLDIR" != "" ]; then
    CONFIGURE_FLAGS="${CONFIGURE_FLAGS} --with-gmp=$GMP_INSTALLDIR"
fi

if [ "x$NO_OPENMP" != x ]; then
    CONFIGURE_FLAGS="${CONFIGURE_FLAGS} --disable-openmp"
fi

LDFLAGS_SET="${LDFLAGS} -L${PREFIX}/lib/"

CPPFLAGS_SET="-I ${PREFIX}/include"

echo "installing shared"

mkdir -p build
cd build

if [ "x$NMZ_EXTENDED_TESTS" != x ]; then
../configure ${CONFIGURE_FLAGS}  CPPFLAGS="-I ${PREFIX}/include -DNMZ_EXTENDED_TESTS" LDFLAGS="${LDFLAGS_SET}"  $EXTRA_FLAGS --srcdir=..
else
../configure ${CONFIGURE_FLAGS}  CPPFLAGS="${CPPFLAGS_SET}" LDFLAGS="${LDFLAGS_SET}" $EXTRA_FLAGS --srcdir=..
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
