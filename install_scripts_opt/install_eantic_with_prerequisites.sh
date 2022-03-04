#!/usr/bin/env bash

set -e

if [ "$OSTYPE" != "msys" ]; then
./install_scripts_opt/install_nmz_mpfr.sh
./install_scripts_opt/install_nmz_flint.sh
fi
./install_scripts_opt/install_nmz_arb.sh
./install_scripts_opt/install_nmz_antic.sh
./install_scripts_opt/install_nmz_e-antic.sh
