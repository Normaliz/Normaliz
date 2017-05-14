#! /bin/sh
PYNORMALIZ_DIR=`pwd`
set -e
NORMALIZ_VERSION=3.3.0
NORMALIZ_PREFIX=`pwd`/local
wget https://github.com/Normaliz/Normaliz/releases/download/v$NORMALIZ_VERSION/Normaliz-$NORMALIZ_VERSION.tar.gz
rm -Rf Normaliz-$NORMALIZ_VERSION
tar xfz Normaliz-$NORMALIZ_VERSION.tar.gz
cd normaliz-$NORMALIZ_VERSION
./configure --prefix=$NORMALIZ_PREFIX
make -j4
# Limit number of threads
OMP_NUM_THREADS=4
export OMP_NUM_THREADS
make -j4 check
make install
# Install PyNormaliz
cd $PYNORMALIZ_DIR
export LDFLAGS="-L$NORMALIZ_PREFIX/lib -Wl,-rpath,$NORMALIZ_PREFIX/lib $LDFLAGS"
export CPATH="$NORMALIZ_PREFIX/include:$CPATH"
pip install .
