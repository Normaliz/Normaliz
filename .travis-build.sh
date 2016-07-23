#! /bin/sh
case $BUILDSYSTEM in
    cmake)
	mkdir -p BUILD || exit 1
	cd BUILD || exit 1
	pwd
	cmake ../source || exit 1
	make -j2 || exit 1
	OMP_NUM_THREADS=4 make check || exit 1
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
	NMZDIR=`pwd`
	rm -Rf CoCoA
	mkdir CoCoA || exit 1
	cd CoCoA || exit 1
	wget http://cocoa.dima.unige.it/cocoalib/tgz/CoCoALib-$COCOALIB_VERSION.tgz || exit 1
	tar xf CoCoALib-$COCOALIB_VERSION.tgz || exit 1
	cd CoCoALib-$COCOALIB_VERSION || exit 1
	COCOALIB_DIR=`pwd` || exit 1
	./configure --threadsafe-hack || exit 1
	make -j2 library || exit 1
	cd $NMZDIR || exit 1
	./bootstrap.sh || exit 1
	# Don't pass CoCoA flags here. We want to make sure that the distribution
	# is complete even when this source tree is not configured with nmzintegrate.
	./configure || exit 1
	# Rather, build the unpacked distribution with CoCoA.
	make -j2 DISTCHECK_CONFIGURE_FLAGS="--with-cocoalib=$COCOALIB_DIR --enable-nmzintegrate" distcheck || ( echo '#### Contents of config.log: ####'; cat _build/config.log; exit 1)
	;;
    *)
	# autotools
	./bootstrap.sh || exit 1
	./configure || ( echo '#### Contents of config.log: ####'; cat config.log; exit 1)
	make -j2 || exit 1
	OMP_NUM_THREADS=4 make -j2 check || exit 1
	;;
esac
