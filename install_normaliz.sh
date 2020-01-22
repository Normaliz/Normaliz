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

echo "installing shared"

mkdir -p build_shared
cd build_shared

CPPFLAGS="${CPPFLAGS} -I${PREFIX}/include"
if [ "x$NMZ_EXTENDED_TESTS" != x ]; then
    CPPFLAGS="${CPPFLAGS} -DNMZ_EXTENDED_TESTS"
fi
LDFLAGS="${LDFLAGS} -L${PREFIX}/lib/"

../configure ${CONFIGURE_FLAGS} $EXTRA_FLAGS

make clean
make -j4
make install

# make distclean

cd ..
# rm -r build_shared

cp -f ${PREFIX}/bin/* .
cp ${PREFIX}/lib/libnormaliz.a source/libnormaliz ## for compatibility with Makefile.classic

echo "Normaliz installation complete"
