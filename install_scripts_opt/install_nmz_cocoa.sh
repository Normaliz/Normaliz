#!/usr/bin/env bash

set -e

echo "::group::cocoa"

source $(dirname "$0")/common.sh

if [ "$GMP_INSTALLDIR" != "" ]; then
    CPPFLAGS="${CPPFLAGS} -I${GMP_INSTALLDIR}/include"
    LDFLAGS="${LDFLAGS} -L${GMP_INSTALLDIR}/lib"
fi

if [[ $OSTYPE == darwin* ]] && [ "$GMP_INSTALLDIR" != "" ]; then
    COCOA_CONFIGURE_FLAGS=" --with-libgmp=${GMP_INSTALLDIR}/lib/libgmp.a"
fi

##  script for the installation of CoCoALib
## as far as needed by libnormaliz

COCOA_VERSION="0.99710"
COCOA_URL="https://cocoa.dima.unige.it/cocoa/cocoalib/tgz/CoCoALib-${COCOA_VERSION}.tgz"
COCOA_SHA256=80d472fd74c7972f8f2a239679e7ad8ae8a43676e3c259c2218ae2480a6267a8

echo "Installing CoCoA..."

# download & extract
mkdir -p ${NMZ_OPT_DIR}/CoCoA_source/
cd ${NMZ_OPT_DIR}/CoCoA_source/
if [ "$OSTYPE" != "msys" ]; then
	../../download.sh ${COCOA_URL} ${COCOA_SHA256}
	if [ ! -d CoCoALib-${COCOA_VERSION} ]; then
		tar xvf CoCoALib-${COCOA_VERSION}.tgz
fi
else # pre release version for MSYS
	COCOA_URL=https://github.com/Normaliz/Normaliz/releases/download/v3.9.2/Prerelease_CoCoA_for.MSYS.tgz
	COCOA_VERSION="0.99719"
	../../download.sh ${COCOA_URL}
	tar xvf Prerelease_CoCoA_for.MSYS.tgz
	cd  CoCoALib-${COCOA_VERSION}
	cp ../../../install_scripts_opt/cocoa_patches/SignalWatcher.C src/AlgebraicCore
	cp ../../../install_scripts_opt/cocoa_patches/SignalWatcher.H include/CoCoA
	cp ../../../install_scripts_opt/cocoa_patches/configure .
	cd ..
fi

# configure & compile
cd CoCoALib-${COCOA_VERSION}
if [ ! -f configuration/autoconf.mk ]; then
    ./configure --threadsafe-hack --no-boost ${COCOA_CONFIGURE_FLAGS}
fi
make library -j4
mkdir -p ${PREFIX}/include/CoCoA
cp include/CoCoA/*.H ${PREFIX}/include/CoCoA
# mkdir -p ${PREFIX}/lib ## in common.sh
cp lib/libcocoa.a ${PREFIX}/lib

echo "CoCoALib installed"
