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
        clang++ --version
        ./install_nmz_flint.sh
	;;
esac
# Set up E-ANTIC and dependencies if necessary.
case $BUILDSYSTEM in
    *-enfnormaliz*)
        export NMZ_COMPILER=$CXX
        export PATH="`brew --prefix`/opt/llvm/bin/:$PATH"
        export LDFLAGS="-L`brew --prefix`/opt/llvm/lib"
        clang++ --version
        ./install_nmz_flint_for_eantic.sh > /dev/null && ./install_nmz_arb.sh > /dev/null && ./install_nmz_antic.sh > /dev/null && ./install_nmz_e-antic.sh
        ;;
esac
# Set up CoCoA if necessary for this build.
case $BUILDSYSTEM in
    *-nmzintegrate*)
        export NMZ_COMPILER=$CXX
        clang++ --version
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
    	
        ./configure $CONFIGURE_FLAGS  --prefix=${INSTALLDIR} --with-cocoalib=${INSTALLDIR} --with-flint=${INSTALLDIR} --disable-shared
        
        mkdir -p ${OPTLIBDIR}/hide
        
        ls ${OPTLIBDIR}
        
        if [-f *.dylib ]; then
                mv -f ${OPTLIBDIR}/*.dylib.* ${OPTLIBDIR}/hide
                mv -f ${OPTLIBDIR}/*.dylib ${OPTLIBDIR}/hide
                mv -f ${OPTLIBDIR}/*la ${OPTLIBDIR}/hide
        fi        
        if [ -f *.so ]; then
                mv -f ${OPTLIBDIR}/*.so.* ${OPTLIBDIR}/hide
                mv -f ${OPTLIBDIR}/*.so ${OPTLIBDIR}/hide
                mv -f ${OPTLIBDIR}/*la ${OPTLIBDIR}/hide
        fi
        
        make -j2
        make install
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
