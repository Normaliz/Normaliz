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

E_ANTIC_VERSION=1.1.0
E_ANTIC_URL="https://github.com/flatsurf/e-antic/releases/download/${E_ANTIC_VERSION}/e-antic-${E_ANTIC_VERSION}.tar.gz"
E_ANTIC_SHA256=b914d082c486a9d17cd46b38a55170ba7b1838a2ab8a29f986899d2668f8210e

CONFIGURE_FLAGS="${CONFIGURE_FLAGS} --prefix=${PREFIX} --disable-silent-rules --without-byexample --without-doc --without-benchmark --without-pyeantic"

echo "Installing E-ANTIC..."

mkdir -p ${NMZ_OPT_DIR}/E-ANTIC_source/
cd ${NMZ_OPT_DIR}/E-ANTIC_source

../../download.sh ${E_ANTIC_URL} ${E_ANTIC_SHA256}
if [ ! -d e-antic-${E_ANTIC_VERSION} ]; then
    tar -xvf e-antic-${E_ANTIC_VERSION}.tar.gz
	cd e-antic-${E_ANTIC_VERSION}/libeantic
	sed -i -e s/fmpq_poly_add_fmpq/fmpq_poly_add_fmpq_eantic/g upstream/patched/fmpq_poly_add_fmpq.c
	sed -i -e s/fmpq_poly_add_fmpq/fmpq_poly_add_fmpq_eantic/g upstream/patched/nf_elem_add_fmpq.c
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
