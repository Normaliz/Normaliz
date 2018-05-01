#!/usr/bin/env bash

set -e

WITH_GMP=""
if [ "$GMP_INSTALLDIR" != "" ]; then
  WITH_GMP="--with-gmp=$GMP_INSTALLDIR"
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

## script for the installation of ANTIC for the use in libnormaliz

E_ANTIC_BRANCH=winfried
E_ANTIC_COMMIT=2bfc74e4244ef8eb33cb00e050ae71fc6358abb5

PREFIX=${NMZ_OPT_DIR}

echo "Installing E-ANTIC..."

mkdir -p ${NMZ_OPT_DIR}/E-ANTIC_source/
cd ${NMZ_OPT_DIR}/E-ANTIC_source
if [ -d e-antic ]; then
    (cd e-antic && git fetch origin ${E_ANTIC_BRANCH})
else
    git clone --branch=${E_ANTIC_BRANCH} --single-branch https://github.com/mkoeppe/e-antic
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
export LDFLAGS="-L${NMZ_OPT_DIR}/lib ${LDFLAGS}"
if [ ! -f configure ]; then
    ./bootstrap.sh
fi
if [ ! -f config.status ]; then
    ./configure --prefix=${PREFIX} $WITH_GMP --with-flint="${NMZ_OPT_DIR}" \
                --with-mpfr="${NMZ_OPT_DIR}"
fi
make -j4
make install
