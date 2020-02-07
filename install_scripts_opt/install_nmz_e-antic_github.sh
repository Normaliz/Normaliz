#!/usr/bin/env bash

set -e

source $(dirname "$0")/common.sh

CONFIGURE_FLAGS="--prefix=${PREFIX}"
if [ "x$NO_OPENMP" != x ]; then
    CONFIGURE_FLAGS="${CONFIGURE_FLAGS} --disable-openmp"
fi

if [ "$GMP_INSTALLDIR" != "" ]; then
    CPPFLAGS="${CPFFLAGS} -I${GMP_INSTALLDIR}/include"
    LDFLAGS="${LDFLAGS} -L${GMP_INSTALLDIR}/lib"
fi

## script for the installation of e-antic for the use in libnormaliz

E_ANTIC_BRANCH=master
E_ANTIC_COMMIT=561fb96cdfede786250dd743eb4e2ece182636b8

if [ "x$NO_OPENMP" != x ]; then
    CONFIGURE_FLAGS="${CONFIGURE_FLAGS} --disable-openmp"
fi

echo "Installing E-ANTIC..."

# download & extract
mkdir -p ${NMZ_OPT_DIR}/E-ANTIC_source/
cd ${NMZ_OPT_DIR}/E-ANTIC_source
if [ -d e-antic ]; then
    (cd e-antic && git fetch origin ${E_ANTIC_BRANCH})
else
    git clone --branch=${E_ANTIC_BRANCH} --single-branch https://github.com/videlec/e-antic
fi
cd e-antic
if [ -n "${E_ANTIC_COMMIT}" ]; then
    git checkout ${E_ANTIC_COMMIT}
else
    git pull --ff-only
fi

if [ ! -f configure ]; then
    ./bootstrap.sh
fi

# configure & compile
if [ ! -f config.status ]; then
    ./configure --prefix=${PREFIX} ${BLOCK_OPENMP} \
              CPPFLAGS="${CPPFLAGS} -fPIC"
# --enable-flint-devel ## for Flint development version
fi
make -j4
make install
