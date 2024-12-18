#!/bin/bash
set -e # exit on errors
set -x # print commands and their arguments as they are executed

if [ -e VARS1 ]; then
  source VARS1
fi

echo "CLANG CLANG"
clang++ --version
echo "G++ G++"
g++ --version

case $BUILDSYSTEM in

    *static*)
        OPTLIBDIR=${PREFIX}/lib

        echo "&&&&&&&&&&&&"
        echo ${PREFIX}
        echo ${OPTLIBDIR}
        echo "&&&&&&&&&&&"

        # Remove shared libraries and libtool *.la files to force static linking
        # ls -laR ${OPTLIBDIR}
        rm -f ${OPTLIBDIR}/*.dylib*
        rm -f ${OPTLIBDIR}/*.so*
        rm -f ${OPTLIBDIR}/*.la
        # ls -laR ${OPTLIBDIR}
        if [[ $OSTYPE == darwin* ]]; then
            BREWDIR=$(brew --prefix)
            rm -f ${BREWDIR}/lib/*gmp*.dylib*
            rm -f ${BREWDIR}/lib/*mpfr*.dylib*
            rm -f ${BREWDIR}/lib/*flint*.dylib*
            rm -f ${BREWDIR}/opt/libomp/lib/*.dylib
            rm -f  ${BREWDIR}/opt/llvm/lib/libomp.dylib
            LDFLAGS="${LDFLAGS} -L$(brew --prefix)/opt/libomp/lib"
        fi

        echo "++++++++++"
        echo ${LDFLAGS}
        echo "++++++++++"

        make -j2 LDFLAGS="${LDFLAGS}"
        make install
        ;;

    *makedistcheck*)
        # make -j2 distcheck # checks in succeeding script
        ;;

    *)
        make -j2 -k
        # make -j2 -k check
        # make install
        # if [[ $OSTYPE == darwin* ]]; then
        #     otool -L ${PREFIX}/bin/*
        # else
        #     ldd ${PREFIX}/bin/*
        # fi
        # make installcheck
        ;;
esac

export -p | sed 's/declare -x/export/g' > VARS2
echo "environment variables stored in VARS2 -- this file gets sourced by the succeeding script"
