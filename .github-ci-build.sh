#!/bin/bash
set -e # exit on errors
set -x # print commands and their arguments as they are executed

if [[ $BUILDSYSTEM != *static* ]]; then
    if [[ $OSTYPE == darwin* ]]; then 
        export NMZ_COMPILER=clang++
    fi
fi

source install_scripts_opt/common.sh

# Prepare configure flags
CONFIGURE_FLAGS="${CONFIGURE_FLAGS} --prefix=${PREFIX}"

if [ "x$NO_OPENMP" != x ]; then
    CONFIGURE_FLAGS="${CONFIGURE_FLAGS} --disable-openmp"
fi

# install dependencies
if [[ $BUILDSYSTEM == *nauty* ]]; then
    CONFIGURE_FLAGS="${CONFIGURE_FLAGS} --with-nauty"
else
    CONFIGURE_FLAGS="${CONFIGURE_FLAGS} --without-nauty"
fi

if [[ $BUILDSYSTEM == *flint* || $BUILDSYSTEM == *eantic* ]]; then
    CONFIGURE_FLAGS="${CONFIGURE_FLAGS} --with-flint"
else
    CONFIGURE_FLAGS="${CONFIGURE_FLAGS} --without-flint"
fi

if [[ $BUILDSYSTEM == *eantic* ]]; then
    CONFIGURE_FLAGS="${CONFIGURE_FLAGS} --with-e-antic"
else
    CONFIGURE_FLAGS="${CONFIGURE_FLAGS} --without-e-antic"
fi

if [[ $BUILDSYSTEM == *cocoa* ]]; then
    CONFIGURE_FLAGS="${CONFIGURE_FLAGS} --with-cocoalib"
else
    CONFIGURE_FLAGS="${CONFIGURE_FLAGS} --without-cocoalib"
fi

# Build Normaliz.
./bootstrap.sh

if [[ -z $NO_COVERAGE ]]; then
    export CFLAGS="${CFLAGS} --coverage -O2 -g"
    export CXXFLAGS="${CXXFLAGS} --coverage -O2 -g"
    export LDFLAGS="${LDFLAGS} --coverage"
fi

if [[ $BUILDSYSTEM == *extended* ]]; then
    export CPPFLAGS="${CPPFLAGS} -DNMZ_EXTENDED_TESTS"
fi

./configure ${CONFIGURE_FLAGS} || ( echo '#### Contents of config.log: ####'; cat config.log; exit 1)

# Have normaliz testsuite print running time:
export NICE=time

# Limit number of threads
export OMP_NUM_THREADS=4

case $BUILDSYSTEM in

    *static*)
        OPTLIBDIR=${PREFIX}/lib

        # Remove shared libraries and libtool *.la files to force static linking
        #rm -f ${OPTLIBDIR}/*.dylib*
        #rm -f ${OPTLIBDIR}/*.so*
        #rm -f ${OPTLIBDIR}/*la
        if [[ $OSTYPE == darwin* ]]; then
            BREWDIR=$(brew --prefix)
            rm -f ${BREWDIR}/lib/*gmp*.dylib*
            rm -f ${BREWDIR}/lib/*mpfr*.dylib*
            rm -f ${BREWDIR}/lib/*flint*.dylib*
            export CPPFLAGS="${CPPFLAGS} -I/Library/Developer/CommandLineTools/SDKs/MacOSX11.0.sdk/usr/include/machine"
            # export CPPFLAGS="${CPPFLAGS} -I~/usr/include/machine"
        fi

        make -j2 LDFLAGS="${LDFLAGS} -all-static"
        make install

        if [[ $OSTYPE == darwin* ]]; then
            install -m 0644 /usr/local/opt/llvm/lib/libomp.dylib ${PREFIX}/bin
            install_name_tool -id "@loader_path/./libomp.dylib" ${PREFIX}/bin/libomp.dylib
            install_name_tool -change "/usr/local/opt/llvm/lib/libomp.dylib" "@loader_path/./libomp.dylib" ${PREFIX}/bin/normaliz
        fi

        if [[ $OSTYPE == darwin* ]]; then
            otool -L ${PREFIX}/bin/*
        else
            ldd ${PREFIX}/bin/*
        fi
        
        export NORMPARA=-x=1 ## more is not possible on GitHub Mac

        make check
        ;;

    *makedistcheck*)
        make -j2 distcheck
        ;;

    *)
        make -j2 -k
        make -j2 -k check
        make install
        if [[ $OSTYPE == darwin* ]]; then
            otool -L ${PREFIX}/bin/*
            export CPPFLAGS="${CPPFLAGS} -I/Library/Developer/CommandLineTools/SDKs/MacOSX11.0.sdk/usr/include/machine"
            # export CPPFLAGS="${CPPFLAGS} -I~/usr/include/machine"
        else
            ldd ${PREFIX}/bin/*
        fi
        make installcheck
        ;;
esac
