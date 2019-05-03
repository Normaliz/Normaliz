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

# TODO: directory tools? playground? linear_algebra_test?
# TODO: jNormaliz? java libs in lib/? documentation (papers)?

PYNORMALIZ_COMMIT_SHA=1d98801c98ea91faf2150dc60cbdda3a64df6640

dist-hook:
	git clone https://github.com/Normaliz/PyNormaliz $(distdir)/PyNormaliz
	cd $(distdir)/PyNormaliz && git checkout $(PYNORMALIZ_COMMIT_SHA)
	rm -rf $(distdir)/PyNormaliz/.git

