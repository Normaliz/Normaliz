#! /bin/sh
### Configure, make, make install to prepare an "almost static" binary distribution for Mac OS X.
### Requires Homebrew.
set -e
mkdir -p lib
install -m 0644 `brew --prefix`/opt/gmp/lib/libgmp*.a lib/
rm -Rf local-bindist
../configure --srcdir=.. --with-cocoalib=/Users/mkoeppe/s/sage/Normaliz/CoCoA/CoCoALib-0.99543 --prefix=`pwd`/local-bindist --disable-shared  CC=clang CXX=clang++ LDFLAGS=-L`pwd`/lib
make -j8
make install
### Next step: Run .make-bindist.sh
