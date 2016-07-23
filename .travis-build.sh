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
    *)
	# autotools
	./bootstrap.sh || exit 1
	./configure || ( echo Contents of config.log:; cat config.log; exit 1)
	make -j2 || exit 1
	OMP_NUM_THREADS=4 make -j2 check || exit 1
	;;
esac
