#!/usr/bin/env bash

set -e

if [ "x$NMZ_OPT_DIR" = x ]; then 
    export NMZ_OPT_DIR="${PWD}"/nmz_opt_lib
        mkdir -p ${NMZ_OPT_DIR}
fi

if [ "x$NMZ_COMPILER" != x ]; then
    export CXX="$NMZ_COMPILER"
elif [[ $OSTYPE == darwin* ]]; then
    export NMZ_COMPILER=clang++
    export PATH="`brew --prefix`/opt/llvm/bin/:$PATH"
    export LDFLAGS="-L`brew --prefix`/opt/llvm/lib"
fi
./install_nmz_cocoa.sh
./install_nmz_flint_for_eantic.sh
./install_nmz_arb.sh
./install_nmz_antic.sh
./install_nmz_e-antic.sh
./install_normaliz.sh
