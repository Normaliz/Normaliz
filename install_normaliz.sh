#!/usr/bin/env bash

set -e

./bootstrap.sh

source $(dirname "$0")/install_scripts_opt/common.sh


CONFIGURE_FLAGS="--prefix=${PREFIX}"
if [ "$GMP_INSTALLDIR" != "" ]; then
    CONFIGURE_FLAGS="${CONFIGURE_FLAGS} --with-gmp=$GMP_INSTALLDIR"
fi

if [ "x$NO_OPENMP" != x ]; then
    CONFIGURE_FLAGS="${CONFIGURE_FLAGS} --disable-openmp"
fi

if [ "x$NMZ_EXTENDED_TESTS" != x ]; then
    CPPFLAGS="${CPPFLAGS} -DNMZ_EXTENDED_TESTS"
fi

echo "installing shared"

mkdir -p build
cd build

../configure ${CONFIGURE_FLAGS} CPPFLAGS="${CPPFLAGS}" LDFLAGS="${LDFLAGS} -Wl,-s" $EXTRA_FLAGS --srcdir=..

make clean
make -j4
make install
if [ "x$NMZ_SHARED" == x ]; then
    rm source/normaliz
    make -j4 LDFLAGS="${LDFLAGS} -all-static"
    make install
    if [[ $OSTYPE != darwin* ]]; then
        strip --strip-unneeded --remove-section=.comment --remove-section=.note ${PREFIX}/lib/libnormaliz.a
    fi
fi
# make distclean

cd ..
# rm -r build

cp -f ${PREFIX}/bin/* .
cp ${PREFIX}/lib/libnormaliz.a source/libnormaliz ## for compatibility with Makefile.classic

echo "Normaliz installation complete"
