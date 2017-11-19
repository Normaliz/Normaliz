## script for the installation of Flint for the use in libnormaliz
## including the installation of MPFR (needed for Flint)

FLINT_VERSION="2.5.2"
MPFR_VERSION="3.1.6"

mkdir -p ~/MPFR/
cd ~/MPFR
wget http://www.mpfr.org/mpfr-current/mpfr-${MPFR_VERSION}.tar.gz
tar -xvf ${MPFR_VERSION}.tar.gz
cd mpfr-${MPFR_VERSION}
./configure
make -j4
make install

mkdir -p ~/Flint/
cd ~/Flint
wget http://www.flintlib.org/flint-${FLINT_VERSION}.tar.gz
tar -xvf flint-${FLINT_VERSION}.tar.gz
cd flint-${FLINT_VERSION}
make -j4
make install