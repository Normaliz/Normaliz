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

# install dependencies
case $BUILDSYSTEM in
    *-nauty*)
        export NMZ_COMPILER=$CXX
        ./install_scripts_opt/install_nmz_nauty.sh
        ;;
esac
case $BUILDSYSTEM in
    *-flint*)
        export NMZ_COMPILER=$CXX
        ./install_scripts_opt/install_nmz_flint.sh
        ;;
esac
# Set up E-ANTIC and dependencies if necessary.
case $BUILDSYSTEM in
    *-enfnormaliz*)
        export NMZ_COMPILER=$CXX
        ./install_scripts_opt/install_nmz_flint.sh > /dev/null
        ./install_scripts_opt/install_nmz_arb.sh > /dev/null
        if [ "${CONFIGURE_FLAGS}" = "--disable-openmp" ]; then
            export NO_OPENMP="yes"
        fi
        ./install_scripts_opt/install_nmz_e-antic.sh
        ;;
esac
# Set up CoCoA if necessary for this build.
case $BUILDSYSTEM in
    *-nmzintegrate*)
        export NMZ_COMPILER=$CXX
        ./install_scripts_opt/install_nmz_cocoa.sh
        ;;
esac

# Build Normaliz.
cd $NMZDIR
./bootstrap.sh
CONFIGURE_FLAGS="${CONFIGURE_FLAGS} --prefix=${INSTALLDIR}"
case $BUILDSYSTEM in

    *-enfnormaliz*)

        if [[ $OSTYPE == darwin* ]]; then
            if [[ $BUILDSYSTEM == *static* ]]; then
                install -m 0644 `brew --prefix`/opt/gmp/lib/libgmp*.a ${OPTLIBDIR}
                # export LDFLAGS=-L${OPTLIBDIR}
            fi
        fi

        CONFIGURE_FLAGS="${CONFIGURE_FLAGS} --with-cocoalib=${INSTALLDIR}"
        CONFIGURE_FLAGS="${CONFIGURE_FLAGS} --with-nauty=${INSTALLDIR}"
        CONFIGURE_FLAGS="${CONFIGURE_FLAGS} --with-flint=${INSTALLDIR}"
        CONFIGURE_FLAGS="${CONFIGURE_FLAGS} --disable-shared"
        ./configure ${CONFIGURE_FLAGS} || ( echo '#### Contents of config.log: ####'; cat config.log; exit 1)

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

    autotools-makedistcheck)
        ./configure ${CONFIGURE_FLAGS} || ( echo '#### Contents of config.log: ####'; cat config.log; exit 1)

        make -j2 distcheck

        ;;

    autotools-*)
        CONFIGURE_FLAGS="${CONFIGURE_FLAGS} --with-cocoalib=$INSTALLDIR"
        CONFIGURE_FLAGS="${CONFIGURE_FLAGS} --with-flint=$INSTALLDIR"
        CONFIGURE_FLAGS="${CONFIGURE_FLAGS} --with-nauty=${INSTALLDIR}"
        CONFIGURE_FLAGS="${CONFIGURE_FLAGS} --disable-shared"
        ./configure ${CONFIGURE_FLAGS} || ( echo '#### Contents of config.log: ####'; cat config.log; exit 1)

        make -j2 -k
        make -j2 -k check
        make install
        make installcheck
        ;;

    *)
        # autotools, no libraries
        CONFIGURE_FLAGS="${CONFIGURE_FLAGS} --disable-flint"
        ./configure ${CONFIGURE_FLAGS} || ( echo '#### Contents of config.log: ####'; cat config.log; exit 1)

        make -j2 -k
        make -j2 -k check
        make install
        make installcheck
        ;;
esac
