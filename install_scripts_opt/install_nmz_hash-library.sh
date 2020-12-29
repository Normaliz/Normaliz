#!/usr/bin/env bash

set -e

source $(dirname "$0")/common.sh

if [ "$GMP_INSTALLDIR" != "" ]; then
    CPPFLAGS="${CPPFLAGS} -I${GMP_INSTALLDIR}/include"
    LDFLAGS="${LDFLAGS} -L${GMP_INSTALLDIR}/lib"
fi

##  script for the installation of hash-library
## as far as needed by libnormaliz
## https://create.stephan-brumme.com/hash-library/
## https://github.com/stbrumme/hash-library

HASHLIBRARY_VERSION="7"
HASHLIBRARY_URL="https://create.stephan-brumme.com/hash-library/hash-library.zip"
HASHLIBRARY_SHA256=ac405a57e6b87c6952c3e42cec5395805bc81af05298e11cb41f89a9a7b4ed5d

echo "Installing Hash-Library..."

# download & extract
mkdir -p ${NMZ_OPT_DIR}/Hash-library_source
cd ${NMZ_OPT_DIR}/Hash-library_source/
../../download.sh ${HASHLIBRARY_URL} ${HASHLIBRARY_SHA256}
if [ ! -d hash-library ]; then
    unzip -d hash-library hash-library.zip
fi

# compile (the SHA-256 part of) the library:
cd hash-library
sed -i 's/endian.h/sys\/types.h/g' sha256.cpp
# dynamic:
g++ -Wno-deprecated -Wall -pedantic -O3 -funroll-loops -fPIC -shared -o libsha256.so sha256.cpp
# static:
g++ -Wno-deprecated -Wall -pedantic -O3 -funroll-loops -fPIC -static -c -o libsha256.o sha256.cpp
ar rc libsha256.a libsha256.o

mkdir -p ${PREFIX}/include/hash-library
cp sha256.h ${PREFIX}/include/hash-library
# mkdir -p ${PREFIX}/lib ## in common.sh
cp libsha256.so libsha256.a ${PREFIX}/lib

echo "Hash-Library installed"
