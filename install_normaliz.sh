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
echo "OPT_LIB_DIR"
echo $OPT_LIB_DIR

if [ "x$NMZ_SHARED" = x ]; then

    echo "installing with --disable-shared"

    ./configure --prefix="${PREFIX}" --with-cocoalib="${PREFIX}" --with-flint="${PREFIX}" $EXTRA_FLAGS $WITH_GMP --disable-shared
    
    mkdir -p ${OPTLIBDIR}/hide    
    
    if [[ $OSTYPE == darwin* ]]; then
        for file in ${OPTLIBDIR}/*.dylib*; do
            if [[ -f $file ]]; then
                FOUND=true
                break
            fi
        done
        echo "Found see next line"
        echo $FOUND
        if [ x$FOUND != x ]; then
            echo "Hiding Mac"
            mv -f ${OPTLIBDIR}/*.dylib* ${OPTLIBDIR}/hide
            mv -f ${OPTLIBDIR}/*la ${OPTLIBDIR}/hide
        fi
    else
        for file in ${OPTLIBDIR}/*.so*; do
            if [[ -f $file ]]; then
                FOUND=true
                break
            fi
        done
        echo "Found see next line"
        echo $FOUND
        if [ x$FOUND != x ]; then
            echo "Hiding Linux"
            mv -f ${OPTLIBDIR}/*.so* ${OPTLIBDIR}/hide
            mv -f ${OPTLIBDIR}/*la ${OPTLIBDIR}/hide
        fi
    fi
    
    make clean
    make -j4
    make install
    
    ## we move so and la back to their proper location
    if [ x$FOUND != x ]; then 
        mv -f ${OPTLIBDIR}/hide/* ${OPTLIBDIR}
    fi
    rmdir ${OPTLIBDIR}/hide
    
    mkdir -p ${PREFIX}/bin/hide ## hide the non-shared built binaries
    mv ${PREFIX}/bin/*no* ${PREFIX}/bin/hide
fi

echo "installing shared"

./configure --prefix="${PREFIX}" --with-cocoalib="${PREFIX}" --with-flint="${PREFIX}" $EXTRA_FLAGS $WITH_GMP
make clean
make -j4
make install


if [ "x$NMZ_SHARED" = x ]; then
    mv -f ${PREFIX}/bin/hide/* ${PREFIX}/bin
    rmdir ${PREFIX}/bin/hide
fi

cp -f ${PREFIX}/bin/* .
cp ${PREFIX}/lib/libnormaliz.a source/libnormaliz ## for compatibility with Makefile.classic
cp ${PREFIX}/lib/libQnormaliz.a Qsource/libQnormaliz

make clean

echo "Normaliz installation complete"

