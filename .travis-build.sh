#! /bin/bash
set -e # exit on errors
# Have normaliz testsuite print running time:
NICE=time
export NICE
# Limit number of threads
OMP_NUM_THREADS=4
export OMP_NUM_THREADS
# Top-level directory.
NMZDIR=`pwd`
NMZ_OPT_DIR=${PWD}/nmz_opt_lib
clang++ --version
WITH=""
case $BUILDSYSTEM in
    *-nauty*)
        export NMZ_COMPILER=$CXX
        ./install_scripts_opt/install_nmz_nauty.sh
        WITH="$WITH --with-nauty"
    ;;
esac
case $BUILDSYSTEM in
    *-flint*)
        export NMZ_COMPILER=$CXX
        ./install_scripts_opt/install_nmz_flint.sh
        WITH="$WITH --with-flint"
    ;;
esac
# Set up E-ANTIC and dependencies if necessary.
case $BUILDSYSTEM in
    *-enfnormaliz*)
        export NMZ_COMPILER=$CXX
        ./install_scripts_opt/install_nmz_flint.sh > /dev/null        
        ./install_scripts_opt/install_nmz_arb.sh > /dev/null
        if [ "$CONFIGURE_FLAGS" = "--disable-openmp" ]; then
            export NO_OPENMP="yes"
        fi
        ./install_scripts_opt/install_nmz_e-antic.sh
        WITH="$WITH --with-eantic"
        ;;
esac
# Set up CoCoA if necessary for this build.
case $BUILDSYSTEM in
    *-nmzintegrate*)

        export  NMZ_COMPILER=$CXX
        ./install_scripts_opt/install_nmz_cocoa.sh
        WITH="$WITH --with-cocoalib"
        ;;
esac
# Return to directory
cd $NMZDIR
# Installation directory.
INSTALLDIR="`pwd`/local"
OPTLIBDIR=${INSTALLDIR}/lib
# Build Normaliz.
case $BUILDSYSTEM in
    *-enfnormaliz*)
        ./bootstrap.sh || exit 1
        echo ${INSTALLDIR}

        if [[ $OSTYPE == darwin* ]]; then
            if [[ $BUILDSYSTEM == *static* ]]; then
                install -m 0644 `brew --prefix`/opt/gmp/lib/libgmp*.a ${OPTLIBDIR}
                # export LDFLAGS=-L${OPTLIBDIR}
            fi
        fi
        
        export LDFLAGS=-L${INSTALLDIR}/lib
        export CPPFLAGS=-I${INSTALLDIR}/include
        ./configure $CONFIGURE_FLAGS  --prefix=${INSTALLDIR} $WITH --disable-shared

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
        ./bootstrap.sh || exit 1
        ./configure $CONFIGURE_FLAGS || exit 1
        
        make -j2 distcheck || exit 1

        ;;
    *)
        # autotools, no libraries
        ./bootstrap.sh || exit 1

        export LDFLAGS=-L${INSTALLDIR}/lib
        export CPPFLAGS=-I${INSTALLDIR}/include
        ./configure $CONFIGURE_FLAGS --prefix="$INSTALLDIR" --without-flint || ( echo '#### Contents of config.log: ####'; cat config.log; exit 1)

        make -j2 -k ## || exit 1
        make -j2 -k check ## || exit 1
        make install        
        make installcheck
        ;;
esac
set +e # no exit on errors
