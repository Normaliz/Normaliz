#!/usr/bin/env bash

set -e

./install_scripts_opt/install_nmz_cocoa.sh
if [ "$OSTYPE" != "msys" ]; then
./install_scripts_opt/install_nmz_mpfr.sh
./install_scripts_opt/install_nmz_flint.sh
fi
./install_scripts_opt/install_nmz_nauty.sh
./install_normaliz.sh
