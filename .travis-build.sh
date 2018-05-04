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
case $BUILDSYSTEM in
    *-flint*)
        FLINT_VERSION="2.5.2"
        MPFR_VERSION="4.0.1"

        PREFIX=${NMZ_OPT_DIR}

        mkdir -p ${NMZ_OPT_DIR}/MPFR_source/
        cd ${NMZ_OPT_DIR}/MPFR_source
        wget -N http://www.mpfr.org/mpfr-${MPFR_VERSION}/mpfr-${MPFR_VERSION}.tar.gz
        tar -xvf mpfr-${MPFR_VERSION}.tar.gz
        cd mpfr-${MPFR_VERSION}
        ./configure --prefix=${PREFIX}
        make -j4
        make install
        export MPFR_DIR=${NMZ_OPT_DIR}

        mkdir -p ${NMZ_OPT_DIR}/Flint_source/
        cd ${NMZ_OPT_DIR}/Flint_source
        wget -4 -N http://www.flintlib.org/flint-${FLINT_VERSION}.tar.gz
        tar -xvf flint-${FLINT_VERSION}.tar.gz
        cd flint-${FLINT_VERSION}
        ./configure --prefix=${PREFIX} --with-mpfr=${PREFIX}
        make -j4
        make install
        export FLINT_DIR=${NMZ_OPT_DIR}
	;;
esac
# Set up CoCoA if necessary for this build.
case $BUILDSYSTEM in
    *-nmzintegrate*)
	COCOALIB_VERSION=0.99562
	#rm -Rf CoCoA
	COCOADIR=CoCoA
	COCOALIB_DIR=`pwd`/$COCOADIR/CoCoALib-$COCOALIB_VERSION
	if test ! -f $COCOADIR/CoCoALib-$COCOALIB_VERSION/lib/libcocoa.a ; then
	    mkdir -p $COCOADIR || exit 1
	    cd $COCOADIR || exit 1
	    wget http://cocoa.dima.unige.it/cocoalib/tgz/CoCoALib-$COCOALIB_VERSION.tgz || exit 1
	    rm -Rf CoCoALib-$COCOALIB_VERSION
	    tar xf CoCoALib-$COCOALIB_VERSION.tgz || exit 1
	    cd $COCOALIB_DIR || exit 1
	    ./configure --threadsafe-hack || exit 1
            # Patch out use of Boost.  Otherwise, on Mac OS Travis builds
            # CoCoA finds Boost and libcocoa.a has dependencies on boost libraries.
            # As a result, our detection of libcocoa fails.
            sed -i.orig 's/HAVE_BOOST=yes/HAVE_BOOST=no/;s/BOOST_LDLIBS=.*/BOOST_LDLIBS=/;s/-DCoCoA_WITH_BOOST//;' configuration/autoconf.mk
	    make -j2 library || exit 1
            COCOA_DIR="$COCOALIB_DIR"
            export COCOA_DIR   # for cmake build
	fi
        ;;
esac
# Return to directory
cd $NMZDIR
# Installation directory.
INSTALLDIR="`pwd`/local"
# Build Normaliz.
case $BUILDSYSTEM in
    cmake*)
	mkdir -p BUILD || exit 1
	cd BUILD || exit 1
	pwd
	cmake -DCMAKE_INSTALL_PREFIX:PATH="$INSTALLDIR" ../source || exit 1
	make -j2 || exit 1
	make check || exit 1
        make install
        # Test that installation works (like autotools "make installcheck")
        "$INSTALLDIR"/bin/normaliz --version
	;;
    classic-flint*)
	cat > source/Makefile.configuration <<EOF
CXX = ${CXX-g++}
CXXFLAGS += -static
CXXFLAGS += -std=c++0x -Wall -pedantic -O3 -funroll-loops -fopenmp
GMPFLAGS = -lgmpxx -lgmp
CXXFLAGS += -DNMZ_FLINT -I ${FLINT_DIR}/include
FLINTFLAGS = -L${FLINT_DIR}/lib -lflint -lmpfr
LINKFLAGS += \$(FLINTFLAGS) \$(GMPFLAGS)
INSTALLDIR = $INSTALLDIR
EOF
	make -j2 -C source -f Makefile.classic
	make -C test -f Makefile.classic
        make -C source -f Makefile.classic install
        # Test that installation works (like autotools "make installcheck")
        "$INSTALLDIR"/bin/normaliz --version
	;;
    classic*)
	cat > source/Makefile.configuration <<EOF
CXX = ${CXX-g++}
CXXFLAGS += -std=c++0x -Wall -pedantic -O3 -funroll-loops  -fopenmp
GMPFLAGS = -lgmpxx -lgmp
LINKFLAGS += \$(GMPFLAGS)
INSTALLDIR = $INSTALLDIR
EOF
	make -j2 -C source -f Makefile.classic
	make -C test -k -f Makefile.classic
        # Test that installation works (like autotools "make installcheck")
        make -C source -f Makefile.classic install
        "$INSTALLDIR"/bin/normaliz --version
        cp source/Makefile.configuration Qsource/
	make -j2 -C Qsource -f Makefile.classic
	make -C Qtest -k -f Makefile.classic
        # Test that installation works (like autotools "make installcheck")
        make -C Qsource -f Makefile.classic install
        "$INSTALLDIR"/bin/Qnormaliz --version
	;;
    autotools-makedistcheck)
	./bootstrap.sh || exit 1
	./configure $CONFIGURE_FLAGS || exit 1
	make -j2 distcheck || exit 1
	;;
    autotools-makedistcheck-nmzintegrate)
	# This makes sure that the distribution contains the nmzIntegrate sources
	# and that the distribution correctly builds it when CoCoALib is installed.
	cd $NMZDIR || exit 1
	./bootstrap.sh || exit 1
	# Don't pass CoCoA flags here. We want to make sure that the distribution
	# is complete even when this source tree is not configured with nmzintegrate.
	./configure $CONFIGURE_FLAGS --disable-nmzintegrate --disable-scip || exit 1
	# Rather, build the unpacked distribution with CoCoA.
	make -j2 DISTCHECK_CONFIGURE_FLAGS="$CONFIGURE_FLAGS --with-cocoalib=$COCOALIB_DIR --enable-nmzintegrate --disable-scip --disable-shared" distcheck || ( echo '#### Contents of config.log: ####'; cat config.log; echo '#### Contents of .../_build/.../config.log: ####'; cat normaliz-*/_build/config.log || cat normaliz-*/_build/sub/config.log; exit 1)
	;;
    autotools-flint*)
	./bootstrap.sh || exit 1
	./configure $CONFIGURE_FLAGS --prefix="$INSTALLDIR" --with-flint=$FLINT_DIR  || ( echo '#### Contents of config.log: ####'; cat config.log; exit 1)
	make -j2 -k || exit 1
	make -j2 -k check || exit 1
        make install
        make installcheck
	;;
    autotools-nmzintegrate)
	./bootstrap.sh || exit 1
	./configure $CONFIGURE_FLAGS --prefix="$INSTALLDIR" --with-cocoalib=$COCOALIB_DIR --enable-nmzintegrate || ( echo '#### Contents of config.log: ####'; cat config.log; exit 1)
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
