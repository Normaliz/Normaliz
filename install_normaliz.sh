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
./configure --with-cocoalib=${NMZ_OPT_DIR} --with-flint=${NMZ_OPT_DIR}
make -j4
echo "******************************************************"
echo -e "\033[0;31mIf you want to install normaliz, run sudo make install\033[0m"
echo "******************************************************"
