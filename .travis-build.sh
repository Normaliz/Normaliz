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
# Set up SCIP if necessary for this build
if test x$SCIPOPTSUITE_VERSION = x ; then
    SCIPOPTSUITE_VERSION=3.2.1
fi
case $BUILDSYSTEM in
    *-scipshared*)
	SCIP_FLAVOR=shared
	SCIP_SHARED=true
	;;
    *-scip*)
	SCIP_FLAVOR=static
	SCIP_SHARED=false
	;;
esac
case $BUILDSYSTEM in
    *-scip*)
	SCIP_BUILD_DIR=build_scip_$SCIP_FLAVOR
	if test ! -f $SCIP_BUILD_DIR/scipoptsuite-$SCIPOPTSUITE_VERSION/.completed_build ; then
	    if test ! -f scipoptsuite-$SCIPOPTSUITE_VERSION.tgz ; then
		wget --post-data 'tname=normaliz-travis-ci-build&license=academic' "http://scip.zib.de/download.php?fname=scipoptsuite-$SCIPOPTSUITE_VERSION.tgz" --output-document=scipoptsuite-$SCIPOPTSUITE_VERSION.tgz
	    fi
	    mkdir -p $SCIP_BUILD_DIR
	    cd $SCIP_BUILD_DIR
	    rm -Rf scipoptsuite-$SCIPOPTSUITE_VERSION
	    tar xf $NMZDIR/scipoptsuite-$SCIPOPTSUITE_VERSION.tgz
	    cd scipoptsuite-$SCIPOPTSUITE_VERSION
	    # Pass our own CXXFLAGS because otherwise SoPlex 2.2.1 will use "-std=c++11",
	    # which the old compilers on Travis CI do not support.
	    make -j2 CXXFLAGS="-std=c++0x -fPIC" VERBOSE=true ZLIB=false GMP=false READLINE=false SHARED=$SCIP_SHARED scipoptlib
	    touch .completed_build
	fi
	;;
esac
# Build Normaliz.
cd $NMZDIR
case $BUILDSYSTEM in
    cmake*)
	mkdir -p BUILD || exit 1
	cd BUILD || exit 1
	pwd
	cmake ../source || exit 1
	make -j2 || exit 1
	make check || exit 1
	;;
    classic-scip*)
	cat > source/Makefile.configuration <<EOF
CXX = ${CXX-g++}
CXXFLAGS += -std=c++0x -Wall -pedantic -O3 -funroll-loops -g -fopenmp
SCIPPATH = $NMZDIR/$SCIP_BUILD_DIR/scipoptsuite-$SCIPOPTSUITE_VERSION
GMPFLAGS = -lgmpxx -lgmp
CXXFLAGS += -DNMZ_SCIP -I \$(SCIPPATH)/scip-$SCIPOPTSUITE_VERSION/src
SCIPFLAGS = -L \$(SCIPPATH)/lib -lscipopt-$SCIPOPTSUITE_VERSION.`uname | tr [A-Z] [a-z]`.`uname -m`.gnu.opt -lreadline -lz
LINKFLAGS += \$(SCIPFLAGS) \$(GMPFLAGS)
EOF
	make -j2 -C source -f Makefile.classic
	make -C test -f Makefile.classic
	;;
    classic*)
	cat > source/Makefile.configuration <<EOF
CXX = ${CXX-g++}
CXXFLAGS += -std=c++0x -Wall -pedantic -O3 -funroll-loops -g -fopenmp
GMPFLAGS = -lgmpxx -lgmp
LINKFLAGS += \$(GMPFLAGS)
EOF
	make -j2 -C source -f Makefile.classic
	make -C test -f Makefile.classic
	;;
    autotools-makedistcheck)
	./bootstrap.sh || exit 1
	./configure || exit 1
	make -j2 distcheck || exit 1
	;;
    autotools-makedistcheck-nmzintegrate)
	# This makes sure that the distribution contains the nmzIntegrate sources
	# and that the distribution correctly builds it when CoCoALib is installed.
	COCOALIB_VERSION=0.99543
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
	    make -j2 library || exit 1
	fi
	cd $NMZDIR || exit 1
	./bootstrap.sh || exit 1
	# Don't pass CoCoA flags here. We want to make sure that the distribution
	# is complete even when this source tree is not configured with nmzintegrate.
	./configure --disable-nmzintegrate --disable-scip || exit 1
	# Rather, build the unpacked distribution with CoCoA.
	make -j2 DISTCHECK_CONFIGURE_FLAGS="--with-cocoalib=$COCOALIB_DIR --enable-nmzintegrate --disable-scip --disable-shared" distcheck || ( echo '#### Contents of config.log: ####'; cat config.log; echo '#### Contents of .../_build/.../config.log: ####'; cat normaliz-*/_build/config.log; cat normaliz-*/_build/sub/config.log; exit 1)
	;;
    autotools-scip*)
	./bootstrap.sh || exit 1
	./configure --enable-scip --with-scipoptsuite-src=$SCIP_BUILD_DIR/scipoptsuite-$SCIPOPTSUITE_VERSION || ( echo '#### Contents of config.log: ####'; cat config.log; exit 1)
	make -j2 || exit 1
	make -j2 check || exit 1
	;;
    *)
	# autotools, no SCIP
	./bootstrap.sh || exit 1
	./configure --disable-scip || ( echo '#### Contents of config.log: ####'; cat config.log; exit 1)
	make -j2 || exit 1
	make -j2 check || exit 1
	;;
esac
