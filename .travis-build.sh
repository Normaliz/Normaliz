#!/bin/bash
set -e # exit on errors
set -x # print commands and their arguments as they are executed

# Have normaliz testsuite print running time:
export NICE=time

# Limit number of threads
export OMP_NUM_THREADS=4

# Record various directory paths
NMZDIR=${PWD}
NMZ_OPT_DIR=${PWD}/nmz_opt_lib
INSTALLDIR=${NMZDIR}/local
OPTLIBDIR=${INSTALLDIR}/lib
export NMZ_COMPILER=$CXX

if [ "x$NO_OPENMP" != x ]; then
    CONFIGURE_FLAGS="${CONFIGURE_FLAGS} --disable-openmp"
fi

# install dependencies
case $BUILDSYSTEM in
    *nauty*)
        ./install_scripts_opt/install_nmz_nauty.sh  &> /dev/null # too much output on travis
        CONFIGURE_FLAGS="${CONFIGURE_FLAGS} --with-nauty=${INSTALLDIR}"
        ;;
esac
case $BUILDSYSTEM in
    *flint* | *eantic*)
        ./install_scripts_opt/install_nmz_flint.sh &> /dev/null # too much output on travis
        CONFIGURE_FLAGS="${CONFIGURE_FLAGS} --with-flint=$INSTALLDIR"
        ;;
esac
# Set up E-ANTIC and dependencies if necessary.
case $BUILDSYSTEM in
    *eantic*)
        ./install_scripts_opt/install_nmz_arb.sh &> /dev/null # too much output on travis
        ./install_scripts_opt/install_nmz_e-antic.sh
        CONFIGURE_FLAGS="${CONFIGURE_FLAGS} --with-eantic=$INSTALLDIR"

        if [[ $OSTYPE == darwin* ]]; then
            if [[ $BUILDSYSTEM == *static* ]]; then
                install -m 0644 `brew --prefix`/opt/gmp/lib/libgmp*.a ${OPTLIBDIR}
                # export LDFLAGS=-L${OPTLIBDIR}
            fi
        fi
        ;;
esac
# Set up CoCoA if necessary for this build.
case $BUILDSYSTEM in
    *cocoa*)
        ./install_scripts_opt/install_nmz_cocoa.sh  &> /dev/null # too much output on travis
        CONFIGURE_FLAGS="${CONFIGURE_FLAGS} --with-cocoalib=$INSTALLDIR"
        ;;
esac

# Build Normaliz.
cd $NMZDIR
./bootstrap.sh

if [[ -z $NO_COVERAGE ]]; then
    export CFLAGS="${CFLAGS} --coverage -O2 -g"
    export CXXFLAGS="${CXXFLAGS} --coverage -O2 -g"
    export LDFLAGS="${LDFLAGS} --coverage"
fi

CONFIGURE_FLAGS="${CONFIGURE_FLAGS} --prefix=${INSTALLDIR}"
case $BUILDSYSTEM in
    makedistcheck)
        ;;
    *)
        CONFIGURE_FLAGS="${CONFIGURE_FLAGS} --disable-shared"
        ;;
esac

case $BUILDSYSTEM in
    *extended*)
        CONFIGURE_FLAGS="${CONFIGURE_FLAGS} CPPFLAGS=-DNMZ_EXTENDED_TESTS"
esac;;

./configure ${CONFIGURE_FLAGS} || ( echo '#### Contents of config.log: ####'; cat config.log; exit 1)

case $BUILDSYSTEM in

    *static*)
        mkdir -p ${OPTLIBDIR}/hide
        if [ -f ${OPTLIBDIR}/libflint.dylib ]; then
                echo "Hiding Mac"
                mv -f ${OPTLIBDIR}/*.dylib.* ${OPTLIBDIR}/hide
                mv -f ${OPTLIBDIR}/*.dylib ${OPTLIBDIR}/hide
                mv -f ${OPTLIBDIR}/*la ${OPTLIBDIR}/hide
        fi
        if [ -f ${OPTLIBDIR}/libflint.so ]; then
                echo "Hiding Linux"
                mv -f ${OPTLIBDIR}/*.so.* ${OPTLIBDIR}/hide
                mv -f ${OPTLIBDIR}/*.so ${OPTLIBDIR}/hide
                mv -f ${OPTLIBDIR}/*la ${OPTLIBDIR}/hide
        fi

        make -j2
        make install

        if [[ $OSTYPE == darwin* ]]; then
            if [[ $BUILDSYSTEM == *static* ]]; then
                    install -m 0644 /usr/local/opt/llvm/lib/libomp.dylib ${INSTALLDIR}/bin
                    install_name_tool -id "@loader_path/./libomp.dylib" ${INSTALLDIR}/bin/libomp.dylib
                    install_name_tool -change "/usr/local/opt/llvm/lib/libomp.dylib" "@loader_path/./libomp.dylib" ${INSTALLDIR}/bin/normaliz
            fi
        fi

        if [[ $OSTYPE == darwin* ]]; then
            otool -L ${INSTALLDIR}/bin/*
        else
            ldd ${INSTALLDIR}/bin/*
        fi

        make check
        ;;

    makedistcheck)
        make -j2 distcheck
        ;;

    *)
        make -j2 -k
        make -j2 -k check
        make install
        make installcheck
        ;;
esac
