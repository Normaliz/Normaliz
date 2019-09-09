ACLOCAL_AMFLAGS = -I m4

AUTOMAKE_OPTIONS = dist-zip

SUBDIRS = source test example

EXTRA_DIST = bootstrap.sh CHANGELOG .travis.yml .travis-build.sh

#EXTRA_DIST += Macaulay2/INSTALL Macaulay2/Normaliz.m2
#EXTRA_DIST += Singular/INSTALL Singular/normaliz.lib
EXTRA_DIST += doc/Normaliz.pdf
EXTRA_DIST += doc/Normaliz.tex
EXTRA_DIST += doc/jNormaliz.jpg
# EXTRA_DIST += README.md

EXTRA_DIST += Dockerfile
EXTRA_DIST += Docker/Dockerfile
EXTRA_DIST += install_normaliz.sh
EXTRA_DIST += install_normaliz_with_opt.sh
EXTRA_DIST += install_normaliz_with_eantic.sh

EXTRA_DIST += install_scripts_opt/install_nmz_cocoa.sh
EXTRA_DIST += install_scripts_opt/install_nmz_flint.sh
EXTRA_DIST += install_scripts_opt/install_nmz_arb.sh
EXTRA_DIST += install_scripts_opt/install_nmz_nauty.sh
EXTRA_DIST += install_scripts_opt/install_nmz_e-antic.sh
EXTRA_DIST += install_scripts_opt/install_eantic_with_prerequisites.sh 

EXTRA_DIST += install_pynormaliz.sh

EXTRA_DIST += COPYING INSTALL

EXTRA_DIST += source/dynamic/dynamic.cpp
EXTRA_DIST += source/dynamic/second.mat
EXTRA_DIST += source/dynamic/third.mat

EXTRA_DIST += source/outerpar/outerpar.cpp
EXTRA_DIST += source/outerpar/small_gens.mat


# TODO: directory tools? playground? linear_algebra_test?
# TODO: jNormaliz? java libs in lib/? documentation (papers)?

PYNORMALIZ_COMMIT_SHA=d3d79f80713b19655b219aca7accf01306836efa

dist-hook:
	git clone https://github.com/Normaliz/PyNormaliz $(distdir)/PyNormaliz
	cd $(distdir)/PyNormaliz && git checkout $(PYNORMALIZ_COMMIT_SHA)
	rm -rf $(distdir)/PyNormaliz/.git

