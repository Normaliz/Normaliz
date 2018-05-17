#! /bin/sh
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
case $BUILDSYSTEM in
    *-flint*)
        export NMZ_COMPILER=$CXX
        ./install_nmz_flint.sh
	;;
esac
# Set up E-ANTIC and dependencies if necessary.
case $BUILDSYSTEM in
    *-enfnormaliz*)
        export NMZ_COMPILER=$CXX
        ./install_nmz_flint_for_eantic.sh > /dev/null 
        ./install_nmz_arb.sh > /dev/null
        ./install_nmz_antic.sh > /dev/null
        ./install_nmz_e-antic.sh
        ;;
esac
# Set up CoCoA if necessary for this build.
case $BUILDSYSTEM in
    *-nmzintegrate*)

        export  NMZ_COMPILER=$CXX
	./install_nmz_cocoa.sh
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
            install -m 0644 `brew --prefix`/opt/gmp/lib/libgmp*.a ${OPTLIBDIR}
            export LDFLAGS=-L${OPTLIBDIR}
        fi
        
        if [[ $OSTYPE == darwin* ]]; then
            echo "COmpiler version"
            clang++ --version
        fi
        
        if [[ $COMPILER_OVERRIDE != homebrew-llvm ]]; then
            $CONFIGURE_FLAGS += --disable-shared
        fi
    	
        ./configure $CONFIGURE_FLAGS  --prefix=${INSTALLDIR} --with-cocoalib=${INSTALLDIR} --with-flint=${INSTALLDIR}
        
        mkdir -p ${OPTLIBDIR}/hide

        if [[ $COMPILER_OVERRIDE != homebrew-llvm ]]; then
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
        fi

        make -j2
        make install
        
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
    autotools-flint*)
	./bootstrap.sh || exit 1
	./configure $CONFIGURE_FLAGS --prefix=$INSTALLDIR --with-flint=$INSTALLDIR  --with-cocoalib=$INSTALLDIR || ( echo '#### Contents of config.log: ####'; cat config.log; exit 1)
	
	make -j2 -k || exit 1
	make -j2 -k check || exit 1
        make install        
        make installcheck
	;;
    autotools-nmzintegrate)
	./bootstrap.sh || exit 1
	./configure $CONFIGURE_FLAGS --prefix=$INSTALLDIR --with-cocoalib=$INSTALLDIR --enable-nmzintegrate || ( echo '#### Contents of config.log: ####'; cat config.log; exit 1)
	
	make -j2 -k || exit 1
	make -j2 -k check || exit 1
        make install        
        make installcheck
	;;
    *)
	# autotools, no Flint
	./bootstrap.sh || exit 1
	./configure $CONFIGURE_FLAGS --prefix="$INSTALLDIR" --disable-flint || ( echo '#### Contents of config.log: ####'; cat config.log; exit 1)

	make -j2 -k || exit 1
	make -j2 -k check || exit 1
        make install        
        make installcheck
	;;
esac
set +e # no exit on errors