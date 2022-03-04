#!/usr/bin/env bash

set -e

if [ "$OSTYPE" != "msys" ]; then
	./install_scripts_opt/install_nmz_cocoa.sh
fi
./install_scripts_opt/install_eantic_with_prerequisites.sh
./install_scripts_opt/install_nmz_nauty.sh
./install_normaliz.sh
