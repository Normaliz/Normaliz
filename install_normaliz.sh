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
OPTLIBDIR=${PREFIX}/lib


if [ "x$NMZSHARED" = x ]; then
    ./configure --prefix="${PREFIX}" --with-cocoalib="${PREFIX}" --with-flint="${PREFIX}" $EXTRA_FLAGS $WITH_GMP --disable-shared
    
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
    
    make clean
    make -j4
    make install
    
    ## we move so and la back to their proper location    
    mv -f ${OPTLIBDIR}/hide/* ${OPTLIBDIR}
    rmdir ${OPTLIBDIR}/hide
    
    mkdir -p ${PREFIX}/bin/hide ## hide the non-shared built binaries
    mv ${PREFIX}/bin/*no* ${PREFIX}/bin/hide
fi

./configure --prefix="${PREFIX}" --with-cocoalib="${PREFIX}" --with-flint="${PREFIX}" $EXTRA_FLAGS $WITH_GMP
make clean
make -j4
make install


## move the non-shared binaries back
if [ "x$NMZSHARED" = x ]; then
    mv ${PREFIX}/bin/hide/*no* ${PREFIX}/bin
    rmdir ${PREFIX}/bin/hide
fi

cp -f local/bin/* .

#echo "******************************************************"
#echo -e "\033[0;31mIf you want to install normaliz, run sudo make install\033[0m"
#echo "******************************************************"
