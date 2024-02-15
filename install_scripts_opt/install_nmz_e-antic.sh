#!/usr/bin/env bash

## script for the installation of e-antic for the use in libnormaliz

set -e

echo "::group::e-antic"

source $(dirname "$0")/common.sh

if [ "$GMP_INSTALLDIR" != "" ]; then
    export CPPFLAGS="${CPPFLAGS} -I${GMP_INSTALLDIR}/include"
    export LDFLAGS="${LDFLAGS} -L${GMP_INSTALLDIR}/lib"
fi

if [ "$OSTYPE" == "msys" ]; then
	echo "Hiding libmpfr.la and libflint.a"
	mkdir -p ${PREFIX}/lib/hide
	mv ${PREFIX}/lib/libmpfr.la ${PREFIX}/lib/hide
	mv ${PREFIX}/lib/libflint.a ${PREFIX}/lib/hide
fi

## E_ANTIC_VERSION=1.2.1
E_ANTIC_VERSION=2.0.2
E_ANTIC_URL="https://github.com/flatsurf/e-antic/releases/download/${E_ANTIC_VERSION}/e-antic-${E_ANTIC_VERSION}.tar.gz"
## E_ANTIC_SHA256=a7bfb92620fd7e42a06efbe89e011abee88f4fbd99bcec34fd8300ae9b1cf543
E_ANTIC_SHA256=8328e6490129dfec7f4aa478ebd54dc07686bd5e5e7f5f30dcf20c0f11b67f60

CONFIGURE_FLAGS="${CONFIGURE_FLAGS} --prefix=${PREFIX} --disable-silent-rules --without-byexample --without-doc --without-benchmark --without-pyeantic"

echo "Installing E-ANTIC..."

mkdir -p ${NMZ_OPT_DIR}/E-ANTIC_source/
cd ${NMZ_OPT_DIR}/E-ANTIC_source

../../download.sh ${E_ANTIC_URL} ${E_ANTIC_SHA256}
if [ ! -d e-antic-${E_ANTIC_VERSION} ]; then
    tar -xvf e-antic-${E_ANTIC_VERSION}.tar.gz
	cd e-antic-${E_ANTIC_VERSION}/libeantic
	sed -i -e s/fmpq_poly_add_fmpq/fmpq_poly_add_fmpq_eantic/g upstream/patched/fmpq_poly_add_fmpq.c
        cp ../../../../install_scripts_opt/e-antic_pataches/nf_elem_add_fmpq.c upstream/patched/
	cd ../..
fi

cd e-antic-${E_ANTIC_VERSION}/libeantic


if [ ! -f config.status ]; then
    ./configure ${CONFIGURE_FLAGS}
fi
make -j4
make install

if [ "$OSTYPE" == "msys" ]; then
	echo "Restoring libmpfr.la and libflint.a"
	cp ${PREFIX}/lib/hide/* ${PREFIX}/lib
fi
