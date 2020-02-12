#!/usr/bin/env bash

# Some common definition used by the various install*sh scripts

if [ "x$NMZ_OPT_DIR" = x ]; then
    export NMZ_OPT_DIR="${PWD}"/nmz_opt_lib
    mkdir -p ${NMZ_OPT_DIR}
fi

if [ "x$NMZ_COMPILER" != x ]; then
    export CXX=$NMZ_COMPILER
elif [[ $OSTYPE == darwin* ]]; then
    export CXX=clang++
    export PATH="`brew --prefix`/opt/llvm/bin/:$PATH"
    if [ "x$NMZ_MAC_STATIC" != x ]; then
        install -m 0644 `brew --prefix`/opt/gmp/lib/libgmp*.a ${OPTLIBDIR}
        export LDFLAGS="-L${OPTLIBDIR}"
    fi
    export LDFLAGS="${LDFLAGS} -L`brew --prefix`/opt/llvm/lib"
fi

if [ "x$NMZ_PREFIX" != x ]; then
    mkdir -p ${NMZ_PREFIX}
    PREFIX=${NMZ_PREFIX}
else
    PREFIX=${PWD}/local
fi
