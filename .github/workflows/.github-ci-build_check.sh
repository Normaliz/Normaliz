#!/bin/bash
set -e # exit on errors
set -x # print commands and their arguments as they are executed

if [ -e VARS2 ]; then
  source VARS2
fi

case $BUILDSYSTEM in

    *static*)
        ## export NORMPARA=-x=1 ## if paralleization should not work

        make -j2 -k check
        make install
        if [[ $OSTYPE == darwin* ]]; then
            echo "WWWWWWWWWWWWWWWWWWWWWW"
            diff --version
            echo "VVVVVVVVVVVVVVVVVVVVVV"
            otool -L ${PREFIX}/bin/*
            echo "UUUUUUUUUUUUUUUUUUUUUU"
        else
            ldd ${PREFIX}/bin/*
        fi

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
