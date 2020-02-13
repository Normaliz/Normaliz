#!/usr/bin/env bash

set -e

./install_scripts_opt/install_nmz_cocoa.sh
./install_scripts_opt/install_nmz_mpfr.sh
./install_scripts_opt/install_nmz_flint.sh
./install_scripts_opt/install_nmz_nauty.sh
./install_normaliz.sh
