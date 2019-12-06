#!/usr/bin/env bash

set -e

WITH_GMP=""
if [ "$GMP_INSTALLDIR" != "" ]; then
    CPPFLAGS="${CPFFLAGS} -I${GMP_INSTALLDIR}/include"
    LDFLAGS="${LDFLAGS} -L${GMP_INSTALLDIR}/lib"
fi

source $(dirname "$0")/common.sh

## script for the installation of e-antic for the use in libnormaliz

E_ANTIC_BRANCH=master
E_ANTIC_COMMIT=561fb96cdfede786250dd743eb4e2ece182636b8

if [ "x$NO_OPENMP" != x ]; then
    export BLOCK_OPENMP="--disable-openmp"
fi

echo "Installing E-ANTIC..."

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
# (In particular on Mac OS X, make sure that our version of MPFR comes
# first in the -L search path, not the one from LLVM or elsewhere.
# E_ANTIC's configure puts it last.)
## export LDFLAGS="-L${NMZ_OPT_DIR}/lib ${LDFLAGS}"
## export LDFLAGS="-L${PREFIX}/lib ${LDFLAGS}"
if [ ! -f configure ]; then
    ./bootstrap.sh
fi
if [ ! -f config.status ]; then
    ./configure --prefix=${PREFIX} ${BLOCK_OPENMP} \
              CFLAGS="${CFLAGS} -I${PREFIX}/include" \
              CPPFLAGS="${CPPFLAGS} -I${PREFIX}/include -fPIC" \
              LDFLAGS="${LDFLAGS} -L/${PREFIX}/lib"
# --enable-flint-devel ## for Flint development version
fi
make -j4
make install
