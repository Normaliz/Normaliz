#!/usr/bin/env bash

set -e


PYTHON2=auto
PYTHON3=auto

while true; do
  case "$1" in
    --python2 ) PYTHON2="$2"; shift 2 ;;
    --python3 ) PYTHON3="$2"; shift 2 ;;
    --prefix ) NMZ_PREFIX="$2"; shift 2 ;;
    --user ) PYTHON_USER="--user"; shift 1 ;;
    --sudo ) SUDO="sudo"; shift 1;;
    -- ) shift; break ;;
    * ) break ;;
  esac
done

if [ "x$NMZ_PREFIX" != x ]; then
    PREFIX=${NMZ_PREFIX}
else
    PREFIX=${PWD}/local
fi

if [ "$PYTHON2" = "auto" ]; then
    PYTHON2=$(command -v python2)
fi

if [ "$PYTHON3" = "auto" ]; then
    PYTHON3=$(command -v python3)
fi

if [ ! -d "PyNormaliz" ]; then
    git clone https://github.com/Normaliz/PyNormaliz
fi

if [ ! -d "PyQNormaliz" ]; then
    git clone https://github.com/Normaliz/PyQNormaliz
fi

cd PyNormaliz

if [ -x "$PYTHON3" ]; then
    $PYTHON3 setup.py build_ext --include-dirs="${PREFIX}/include" --library-dirs="${PREFIX}/lib"
    $SUDO $PYTHON3 setup.py install $PYTHON_USER
fi

if [ -x "$PYTHON2" ]; then
    $PYTHON2 setup.py build_ext --include-dirs="${PREFIX}/include" --library-dirs="${PREFIX}/lib"
    $SUDO $PYTHON2 setup.py install $PYTHON_USER
fi

cd ../PyQNormaliz

if [ -x "$PYTHON3" ]; then
    $PYTHON3 setup.py build_ext --include-dirs="${PREFIX}/include" --library-dirs="${PREFIX}/lib"
    $SUDO $PYTHON3 setup.py install $PYTHON_USER
fi

if [ -x "$PYTHON2" ]; then
    $PYTHON2 setup.py build_ext --include-dirs="${PREFIX}/include" --library-dirs="${PREFIX}/lib"
    $SUDO $PYTHON2 setup.py install $PYTHON_USER
fi
