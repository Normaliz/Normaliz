ACLOCAL_AMFLAGS = -I m4

AUTOMAKE_OPTIONS = dist-zip

SUBDIRS = source test example
SUBDIRS += Qsource Qtest Qexample

EXTRA_DIST = bootstrap.sh CHANGELOG .travis.yml .travis-build.sh

#EXTRA_DIST += Macaulay2/INSTALL Macaulay2/Normaliz.m2
#EXTRA_DIST += Singular/INSTALL Singular/normaliz.lib
EXTRA_DIST += doc/Normaliz.pdf
EXTRA_DIST += doc/Normaliz.tex
EXTRA_DIST += doc/jNormaliz.jpg
# EXTRA_DIST += README.md

EXTRA_DIST += Dockerfile
EXTRA_DIST += Docker/Dockerfile
EXTRA_DIST += install_normaliz_with_opt.sh
EXTRA_DIST += install_nmz_cocoa.sh
EXTRA_DIST += install_nmz_flint.sh
EXTRA_DIST += install_normaliz.sh

EXTRA_DIST += install_nmz_antic.sh
EXTRA_DIST += install_nmz_arb.sh
EXTRA_DIST += install_nmz_e-antic.sh
EXTRA_DIST += install_nmz_flint_for_eantic.sh
EXTRA_DIST += install_normaliz_with_qnormaliz_eantic.sh

EXTRA_DIST += install_pynormaliz.sh

EXTRA_DIST += COPYING INSTALL

# TODO: directory tools? playground? linear_algebra_test?
# TODO: jNormaliz? java libs in lib/? documentation (papers)?

PYNORMALIZ_COMMIT_SHA=eb594e2047c41fa89a5115695220ebcd1f75b822
PYQNORMALIZ_COMMIT_SHA=deac7f345c78fe8dd640e0d5487b52932708e64f

dist-hook:
	git clone https://github.com/Normaliz/PyNormaliz $(distdir)/PyNormaliz
	cd $(distdir)/PyNormaliz && git checkout $(PYNORMALIZ_COMMIT_SHA)
	rm -rf $(distdir)/PyNormaliz/.git
	git clone https://github.com/Normaliz/PyQNormaliz $(distdir)/PyQNormaliz
	cd $(distdir)/PyQNormaliz && git checkout $(PYQNORMALIZ_COMMIT_SHA)
	rm -rf $(distdir)/PyQNormaliz/.git
