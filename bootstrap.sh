#! /bin/sh
## Set up the autotools build system.
## This script is only necessary when Normaliz is checked out from git.
## A distribution made using "make distcheck" will contain the complete
## build system already, so it is ready for "./configure && make && make install".
autoreconf -i --force
echo "Now you can do:"
echo "  ./configure && make && make install"
