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

PREFIX=${PWD}/local

## we hide the shared libraries to make libnormaliz and libQnormaliz independent of them
## by forcing the linker to take *.a
OPTLIBDIR=${PREFIX}/lib
mkdir -p ${OPTLIBDIR}/hide
mv -f ${OPTLIBDIR}/*.so.* ${OPTLIBDIR}/hide
-f ${OPTLIBDIR}/*.so ${OPTLIBDIR}/hide
mv -f ${OPTLIBDIR}/*la ${OPTLIBDIR}/hide


./configure --prefix="${PREFIX}" --with-cocoalib="${PREFIX}" --with-flint="${PREFIX}" $EXTRA_FLAGS $WITH_GMP --disable-shared
make clean
make -j4
make install

## we move so and la back to their proper location

mv -f ${OPTLIBDIR}/hide/*.so.* ${OPTLIBDIR}
mv -f ${OPTLIBDIR}/hide/*.so ${OPTLIBDIR}
mv -f ${OPTLIBDIR}/hide/*la ${OPTLIBDIR}

#echo "******************************************************"
#echo -e "\033[0;31mIf you want to install normaliz, run sudo make install\033[0m"
#echo "******************************************************"
