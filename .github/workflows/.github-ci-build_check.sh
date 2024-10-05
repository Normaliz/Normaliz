#!/bin/bash
set -e # exit on errors
set -x # print commands and their arguments as they are executed

if [ -e VARS2 ]; then
  source VARS2
fi

case $BUILDSYSTEM in

    *static*)
        if [[ $OSTYPE == darwin* ]]; then
            install -m 0644 $(brew --prefix)/opt/libomp/lib/libomp.dylib ${PREFIX}/bin
            otool -L ${PREFIX}/bin/normaliz
            install_name_tool -id "@loader_path/./libomp.dylib" ${PREFIX}/bin/libomp.dylib
            install_name_tool -change "$(brew --prefix)/opt/libomp/lib/libomp.dylib" "@loader_path/./libomp.dylib" ${PREFIX}/bin/normaliz
           otool -L ${PREFIX}/bin/normaliz
         fi

        if [[ $OSTYPE == darwin* ]]; then
            otool -L ${PREFIX}/bin/*
        else
            ldd ${PREFIX}/bin/*
        fi

        ## export NORMPARA=-x=1 ## paralleization does not work at present

        make check

        cd local/bin
        echo "CONTENTS OF LOCAL/BIN"
        ls .
        if [ "x$INTEL" != x ]; then
            zip MacOSbinary-Intel.zip *
        else
            zip MacOSbinary.zip *
        fi

        ;;

    *makedistcheck*)
        make -j2 distcheck
        ;;

    *)
        # make -j2 -k
        make -j2 -k check
        make install
        if [[ $OSTYPE == darwin* ]]; then
            echo "VVVVVVVVVVVVVVVVVVVVVV"
            diff --version
            echo "VVVVVVVVVVVVVVVVVVVVVV"
            otool -L ${PREFIX}/bin/*
        else
            ldd ${PREFIX}/bin/*
        fi
        make installcheck
        ;;
esac

export -p | sed 's/declare -x/export/g' > VARS3
