#!/usr/bin/env bash

set -e

source $(dirname "$0")/install_scripts_opt/common.sh

OPTLIBDIR=${PREFIX}/lib
echo "OPT_LIB_DIR"
echo $OPT_LIB_DIR

WITH_GMP=""
if [ "$GMP_INSTALLDIR" != "" ]; then
  WITH_GMP="--with-gmp=$GMP_INSTALLDIR"
fi

if [ "x$NMZ_COMPILER" != x ]; then
    export CXX="$NMZ_COMPILER"
elif [[ $OSTYPE == darwin* ]]; then
	export CXX=clang++
	export PATH="`brew --prefix`/opt/llvm/bin/:$PATH"
	if [ "x$NMZ_MAC_STATIC" != x ]; then
		install -m 0644 `brew --prefix`/opt/gmp/lib/libgmp*.a ${OPTLIBDIR}
    	export LDFLAGS="-L${OPTLIBDIR} -L`brew --prefix`/opt/llvm/lib"
	else
		export LDFLAGS="-L`brew --prefix`/opt/llvm/lib"
	fi
fi

#mkdir -p BUILD
#cd BUILD
if [ ! -e configure ];
then
    ./bootstrap.sh
fi

if [ "x$NO_OPENMP" != x ]; then
    export BLOCK_OPENMP="--disable-openmp"
fi

if [ "x$NMZ_SHARED" = x ]; then

    echo "installing with --disable-shared"
    
    mkdir -p ${OPTLIBDIR}/save_dynamic
    
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
            cp -f ${OPTLIBDIR}/*.la* ${OPTLIBDIR}/save_dynamic
            cp -f ${OPTLIBDIR}/*.dylib* ${OPTLIBDIR}/save_dynamic
            rm -f ${OPTLIBDIR}/*.dylib*
            rm -f ${OPTLIBDIR}/*la
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
            cp -f ${OPTLIBDIR}/*.so* ${OPTLIBDIR}/save_dynamic
            cp -f ${OPTLIBDIR}/*.la* ${OPTLIBDIR}/save_dynamic
            rm -f ${OPTLIBDIR}/*.so*
            rm -f ${OPTLIBDIR}/*la
        fi
    fi

    
    mkdir -p build_static
    
    cd build_static
    
    ../configure --prefix="${PREFIX}" --with-cocoalib="${PREFIX}" --with-flint="${PREFIX}" $EXTRA_FLAGS $WITH_GMP --disable-shared ${BLOCK_OPENMP} --srcdir=..

    make clean
    make -j4
    make install
    
    #make distclean
    
    cd ..
    #rm -r build_static
    
    ## we copy so and la back to their proper location

    cp -f ${OPTLIBDIR}/save_dynamic/* ${OPTLIBDIR}
    
    mkdir -p ${PREFIX}/bin/hide ## hide the non-shared built binaries
    mv -f ${PREFIX}/bin/*no* ${PREFIX}/bin/hide
fi

echo "installing shared"

mkdir -p build_shared
cd build_shared

../configure --prefix="${PREFIX}" --with-cocoalib="${PREFIX}" --with-flint="${PREFIX}" $EXTRA_FLAGS $WITH_GMP ${BLOCK_OPENMP} --srcdir=..
make clean
make -j4
make install

# make distclean

cd ..
# rm -r build_shared


if [ "x$NMZ_SHARED" = x ]; then
    mv -f ${PREFIX}/bin/hide/* ${PREFIX}/bin
    rmdir ${PREFIX}/bin/hide
fi

if [[ $OSTYPE == darwin* ]]; then
        if [ "x$NMZ_MAC_STATIC" != x ]; then
                install -m 0644 /usr/local/opt/llvm/lib/libomp.dylib ${PREFIX}/bin
                install_name_tool -id "@loader_path/./libomp.dylib" ${PREFIX}/bin/libomp.dylib
                install_name_tool -change "/usr/local/opt/llvm/lib/libomp.dylib" "@loader_path/./libomp.dylib" ${PREFIX}/bin/normaliz
	fi
fi

cp -f ${PREFIX}/bin/* .
cp ${PREFIX}/lib/libnormaliz.a source/libnormaliz ## for compatibility with Makefile.classic

echo "Normaliz installation complete"

