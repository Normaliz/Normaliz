#!/bin/bash
set -e # exit on errors
set -x # print commands and their arguments as they are executed

if [[ $BUILDSYSTEM != *static* ]]; then
    if [[ $OSTYPE == darwin* ]]; then 
        export NMZ_COMPILER=clang++
    fi
fi

source install_scripts_opt/common.sh

# Prepare configure flags
CONFIGURE_FLAGS="${CONFIGURE_FLAGS} --prefix=${PREFIX}"

if [ "x$NO_OPENMP" != x ]; then
    CONFIGURE_FLAGS="${CONFIGURE_FLAGS} --disable-openmp"
fi

# install dependencies
if [[ $BUILDSYSTEM == *nauty* ]]; then
    CONFIGURE_FLAGS="${CONFIGURE_FLAGS} --with-nauty"
else
    CONFIGURE_FLAGS="${CONFIGURE_FLAGS} --without-nauty"
fi

if [[ $BUILDSYSTEM == *flint* || $BUILDSYSTEM == *eantic* ]]; then
    CONFIGURE_FLAGS="${CONFIGURE_FLAGS} --with-flint"
else
    CONFIGURE_FLAGS="${CONFIGURE_FLAGS} --without-flint"
fi

if [[ $BUILDSYSTEM == *eantic* ]]; then
    CONFIGURE_FLAGS="${CONFIGURE_FLAGS} --with-e-antic"
else
    CONFIGURE_FLAGS="${CONFIGURE_FLAGS} --without-e-antic"
fi

if [[ $BUILDSYSTEM == *cocoa* ]]; then
    CONFIGURE_FLAGS="${CONFIGURE_FLAGS} --with-cocoalib"
else
    CONFIGURE_FLAGS="${CONFIGURE_FLAGS} --without-cocoalib"
fi

if [[ $BUILDSYSTEM == *static* ]]; then
    CONFIGURE_FLAGS="${CONFIGURE_FLAGS} --disable-shared"
    export LDFLAGS="${LDFLAGS} -all-static"

    if [[ $BUILDSYSTEM == *eantic* ]]; then
      sed -ie s/-leanticxx\ -leantic/-leanticxx\ -leantic\ -lantic\ -larb/g configure.ac
    fi
fi


# Build Normaliz.
./bootstrap.sh

if [[ -z $NO_COVERAGE ]]; then
    export CFLAGS="${CFLAGS} --coverage -O2 -g"
    export CXXFLAGS="${CXXFLAGS} --coverage -O2 -g"
    export LDFLAGS="${LDFLAGS} --coverage"
fi

if [[ $BUILDSYSTEM == *extended* ]]; then
    export CPPFLAGS="${CPPFLAGS} -DNMZ_EXTENDED_TESTS"
fi

export -p | sed 's/declare -x/export/g' > VARS0
echo "environment variables stored in VARS0 (for debugging)"

./configure ${CONFIGURE_FLAGS} || ( echo '#### Contents of config.log: ####'; cat config.log; exit 1)

# Have normaliz testsuite print running time:
export NICE=time

# Limit number of threads
export OMP_NUM_THREADS=4

export -p | sed 's/declare -x/export/g' > VARS1
echo "environment variables stored in VARS1 -- this file gets sourced by the succeeding script"
