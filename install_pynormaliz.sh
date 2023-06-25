#!/usr/bin/env bash

set -e  # exit script on errors

PYTHON3=auto

while true; do
  case "$1" in
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

if [ "$PYTHON3" = "auto" ]; then
    PYTHON3=$(command -v python3)
fi

if [ ! -d "PyNormaliz" ]; then
    git clone https://github.com/Normaliz/PyNormaliz
fi

cd PyNormaliz


export NORMALIZ_LOCAL_DIR="${PREFIX}"

if [ -x "$PYTHON3" ]; then
    rm -rf build
    # Do not invoke setup.py directly. This is deprecated functionality
    # of setuptools.
    $SUDO $PYTHON3 -m pip install -v -v $PYTHON_USER .
fi
