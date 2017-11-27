#!/bin/sh -e

if [ "x$NMZ_OPT_DIR" = x ]; then 
    export NMZ_OPT_DIR=${PWD}/nmz_opt_lib
        mkdir -p ${NMZ_OPT_DIR}
fi

#mkdir -p BUILD
#cd BUILD
if [ ! -e configure ];
then
    ./bootstrap.sh
fi
export CPPFLAGS="-I ${NMZ_OPT_DIR}/include"
export LDFLAGS="-L${NMZ_OPT_DIR}/lib -Wl,-rpath,${NMZ_OPT_DIR}/lib"
export LIBS="-lflint"
export CXXFLAGS="-DNMZ_FLINT"
export LT_SYS_LIBRARY_PATH="${NMZ_OPT_DIR}/lib"
./configure --with-cocoalib=${NMZ_OPT_DIR}
make -j4
echo "If you want to install normaliz, run sudo make install"