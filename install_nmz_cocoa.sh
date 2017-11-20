##  script for the installation of CoCoALib
## as far as needed by libnormaliz

COCOA_VERSION="0.99550"

INSTALLDIR="/usr/local"

mkdir -p ~/CoCoA/
cd ~/CoCoA/
wget http://cocoa.dima.unige.it/cocoalib/tgz/CoCoALib-${COCOA_VERSION}.tgz
tar xvf CoCoALib-${COCOA_VERSION}.tgz
cd CoCoALib-${COCOA_VERSION}
./configure --threadsafe-hack --no-boost
make library -j4
sudo mkdir -p  ${INSTALLDIR}/include
sudo mkdir  -p ${INSTALLDIR}/include/CoCoA
sudo cp include/CoCoA/*.H  ${INSTALLDIR}/include/CoCoA
sudo mkdir  -p ${INSTALLDIR}/lib
sudo cp lib/libcocoa.a ${INSTALLDIR}/lib