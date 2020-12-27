#!/usr/bin/env bash

# Some common definition used by the various install*sh scripts

if [ "x$NMZ_OPT_DIR" = x ]; then
    export NMZ_OPT_DIR="${PWD}"/nmz_opt_lib
    mkdir -p ${NMZ_OPT_DIR}
fi

if [ "x$NMZ_COMPILER" != x ]; then
    export CXX=$NMZ_COMPILER
elif [[ $OSTYPE == darwin* ]]; then   ## activate Homebrew LLVM
    LLVMDIR="$(brew --prefix)/opt/llvm"
    export LDFLAGS="${LDFLAGS} -L${LLVMDIR}/lib -Wl,-rpath,${LLVMDIR}/lib"
    export CPPFLAGS="${CPPFLAGS} -I ${LLVMDIR}/include"
    export PATH="${LLVMDIR}/bin/:$PATH"
    export CXX=clang++
    echo "CLANG VERSION"
    clang++ --version
fi

if [ "x$NMZ_PREFIX" != x ]; then
    mkdir -p ${NMZ_PREFIX}
    export PREFIX=${NMZ_PREFIX}
else
    export PREFIX=${PWD}/local
fi

if [ "$GMP_INSTALLDIR" != "" ]; then
    export CPPFLAGS="${CPPFLAGS} -I${GMP_INSTALLDIR}/include"
    export LDFLAGS="${LDFLAGS} -L${GMP_INSTALLDIR}/lib"
fi

# Make sure our library versions come first in the search path
export CPPFLAGS="-I${PREFIX}/include ${CPPFLAGS}"
export LDFLAGS="-L${PREFIX}/lib ${LDFLAGS}"

mkdir -p ${PREFIX}/lib
mkdir -p ${PREFIX}/include

echo "**************"
echo $NMZ_OPT_DIR
echo $PREFIX
echo $CPPFLAGS
echo $LDFLAGS
echo "-----------"

