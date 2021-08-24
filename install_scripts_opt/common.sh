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
    export CC=clang
    echo "CLANG++ VERSION"
    clang++ --version
    echo "CLANG VERSION"
    clang --version
fi

if [ "x$NMZ_PREFIX" != x ]; then
    mkdir -p ${NMZ_PREFIX}
    export PREFIX=${NMZ_PREFIX}
else
    export PREFIX=${PWD}/local
fi

if [[ $OSTYPE == darwin* ]] &&  [ "$GMP_INSTALLDIR" == "" ]; then
    GMP_INSTALLDIR=/usr/local
fi

if [ "$GMP_INSTALLDIR" != "" ]; then
    export CPPFLAGS="${CPPFLAGS} -I${GMP_INSTALLDIR}/include"
    export LDFLAGS="${LDFLAGS} -L${GMP_INSTALLDIR}/lib"
fi

# Make sure our library versions come first in the search path
export CPPFLAGS="-I${PREFIX}/include ${CPPFLAGS}"
export LDFLAGS="-L${PREFIX}/lib ${LDFLAGS}"


if [[ $OSTYPE != darwin* ]]; then
    # Since we're installing into a non-standard prefix, we have to help
    # the linker find indirect dependencies such as libantic.so which is a
    # dependency of libeantic.so. (We could also overlink, and link with
    # -lantic but we do not depend on antic directly, so we should not do
    # that; see e.g. http://www.kaizou.org/2015/01/linux-libraries.html.)
    # For some odd reason Debian does not render rpath-link as a RUNPATH in a
    # shared C library, so we set the rpath instead which appears to have the
    # same effect.
    export LDFLAGS="${LDFLAGS} -Wl,-enable-new-dtags -Wl,-rpath=${PREFIX}/lib"
fi

mkdir -p ${PREFIX}/lib
mkdir -p ${PREFIX}/include

echo "**************"
echo $NMZ_OPT_DIR
echo $PREFIX
echo $CPPFLAGS
echo $LDFLAGS
echo "-----------"
