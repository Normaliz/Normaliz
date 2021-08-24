#!/usr/bin/env bash

set -e

echo "::group::e-antic"

source $(dirname "$0")/common.sh

if [ "$GMP_INSTALLDIR" != "" ]; then
    export CPPFLAGS="${CPPFLAGS} -I${GMP_INSTALLDIR}/include"
    export LDFLAGS="${LDFLAGS} -L${GMP_INSTALLDIR}/lib"
fi

## script for the installation of e-antic for the use in libnormaliz
E_ANTIC_VERSION=1.0.3
E_ANTIC_URL="https://github.com/flatsurf/e-antic/releases/download/${E_ANTIC_VERSION}/e-antic-${E_ANTIC_VERSION}.tar.gz"
E_ANTIC_SHA256=eea1dc66fed5962425bc7d2c5ccecb50d25c082b1d84276fa3838bfa96d9cb62

CONFIGURE_FLAGS="--prefix=${PREFIX}"

CONFIGURE_FLAGS="--prefix=${PREFIX} --without-benchmark --disable-silent-rules"

echo "Installing E-ANTIC..."

mkdir -p ${NMZ_OPT_DIR}/E-ANTIC_source/
cd ${NMZ_OPT_DIR}/E-ANTIC_source

../../download.sh ${E_ANTIC_URL} ${E_ANTIC_SHA256}
if [ ! -d e-antic-${E_ANTIC_VERSION} ]; then
    tar -xvf e-antic-${E_ANTIC_VERSION}.tar.gz
fi

cd e-antic-${E_ANTIC_VERSION}/libeantic

sed -i -e s/fmpq_poly_add_fmpq/fmpq_poly_add_fmpq_eantic/g upstream/patched/fmpq_poly_add_fmpq.c
sed -i -e s/fmpq_poly_add_fmpq/fmpq_poly_add_fmpq_eantic/g upstream/patched/nf_elem_add_fmpq.c

if [ ! -f config.status ]; then
    ./configure ${CONFIGURE_FLAGS} --without-byexample --without-doc --without-benchmark --without-pyeantic
fi
make -j4
make install
