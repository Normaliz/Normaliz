#!/usr/bin/env bash

set -e

if [ "x$NMZ_PREFIX" != x ]; then
    mkdir -p ${NMZ_PREFIX}
    PREFIX=${NMZ_PREFIX}
else
    PREFIX=${PWD}/local
fi

OPTLIBDIR=${PREFIX}/lib
echo "OPT_LIB_DIR"
echo $OPT_LIB_DIR

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

