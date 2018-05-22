#!/usr/bin/env bash

set -e

WITH_GMP=""
if [ "$GMP_INSTALLDIR" != "" ]; then
  WITH_GMP="--with-gmp=$GMP_INSTALLDIR"
fi

if [ "x$NMZ_OPT_DIR" = x ]; then 
    export NMZ_OPT_DIR="${PWD}"/nmz_opt_lib
        mkdir -p ${NMZ_OPT_DIR}
fi

if [ "x$NMZ_COMPILER" != x ]; then
    export CXX="$NMZ_COMPILER"
elif [[ $OSTYPE == darwin* ]]; then
    export CXX=clang++
    export PATH="`brew --prefix`/opt/llvm/bin/:$PATH"
    export LDFLAGS="-L`brew --prefix`/opt/llvm/lib"
fi

#mkdir -p BUILD
#cd BUILD
if [ ! -e configure ];
then
    ./bootstrap.sh
fi

if [ "x$NMZ_PREFIX" != x ]; then
    mkdir -p ${NMZ_PREFIX}
    PREFIX=${NMZ_PREFIX}
else
    PREFIX=${PWD}/local
fi
OPTLIBDIR=${PREFIX}/lib

if [ "x$NMZ_SHARED" = x ]; then
    ./configure --prefix="${PREFIX}" --with-cocoalib="${PREFIX}" --with-flint="${PREFIX}" $EXTRA_FLAGS $WITH_GMP --disable-shared
else
    ./configure --prefix="${PREFIX}" --with-cocoalib="${PREFIX}" --with-flint="${PREFIX}" $EXTRA_FLAGS $WITH_GMP
fi

## we hide the shared libraries to make libnormaliz and libQnormaliz independent of them
## by forcing the linker to take *.a

if [ "x$NMZ_SHARED" = x ]; then
    mkdir -p ${OPTLIBDIR}/hide
    if [[ $OSTYPE == darwin* ]]; then
        mv -f ${OPTLIBDIR}/*.dylib.* ${OPTLIBDIR}/hide
        mv -f ${OPTLIBDIR}/*.dylib ${OPTLIBDIR}/hide
        mv -f ${OPTLIBDIR}/*la ${OPTLIBDIR}/hide
    else
        mv -f ${OPTLIBDIR}/*.so.* ${OPTLIBDIR}/hide
        mv -f ${OPTLIBDIR}/*.so ${OPTLIBDIR}/hide
        mv -f ${OPTLIBDIR}/*la ${OPTLIBDIR}/hide
    fi
fi

make clean
make -j4
make install

## we move so and la back to their proper location

if [ "x$NMZ_SHARED" = x ]; then
    mv -f ${OPTLIBDIR}/hide/* ${OPTLIBDIR}
    rmdir ${OPTLIBDIR}/hide
fi

cp -f local/bin/* .

#echo "******************************************************"
#echo -e "\033[0;31mIf you want to install normaliz, run sudo make install\033[0m"
#echo "******************************************************"
