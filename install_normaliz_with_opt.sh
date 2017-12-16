#!/bin/sh -e
if [ "x$NMZ_OPT_DIR" = x ]; then 
    export NMZ_OPT_DIR=${PWD}/nmz_opt_lib
        mkdir -p ${NMZ_OPT_DIR}
fi
./install_nmz_cocoa.sh
./install_nmz_flint.sh
./install_normaliz.sh