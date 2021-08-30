ACLOCAL_AMFLAGS = -I m4

SUBDIRS = source test example

EXTRA_DIST = bootstrap.sh CHANGELOG \
             .github/workflows/ci.yml \
             .github/workflows/.github-ci-build_configure.sh \
             .github/workflows/.github-ci-build_configure.sh \
             .github/workflows/.github-ci-build_configure.sh

EXTRA_DIST += Macaulay2/INSTALL Macaulay2/Normaliz.m2
EXTRA_DIST += Singular/INSTALL Singular/normaliz.lib Singular/normaliz.pdf
EXTRA_DIST += doc/Normaliz.pdf
EXTRA_DIST += doc/Normaliz.tex
EXTRA_DIST += doc/jNormaliz.jpg
# EXTRA_DIST += README.md

EXTRA_DIST += Dockerfile
EXTRA_DIST += Docker/Dockerfile
EXTRA_DIST += download.sh
EXTRA_DIST += install_normaliz.sh
EXTRA_DIST += install_normaliz_with_opt.sh
EXTRA_DIST += install_normaliz_with_eantic.sh
EXTRA_DIST += download.sh

EXTRA_DIST += install_scripts_opt/common.sh
EXTRA_DIST += install_scripts_opt/install_nmz_cocoa.sh
# EXTRA_DIST += install_scripts_opt/install_nmz_cocoa_old.sh
EXTRA_DIST += install_scripts_opt/install_nmz_flint.sh
EXTRA_DIST += install_scripts_opt/install_nmz_mpfr.sh
EXTRA_DIST += install_scripts_opt/install_nmz_arb.sh
EXTRA_DIST += install_scripts_opt/install_nmz_hash-library.sh
EXTRA_DIST += install_scripts_opt/install_nmz_nauty.sh
EXTRA_DIST += install_scripts_opt/install_nmz_e-antic.sh
EXTRA_DIST += install_scripts_opt/install_eantic_with_prerequisites.sh
EXTRA_DIST += install_scripts_opt/install_nmz_antic.sh

EXTRA_DIST += source/chunk/run_single.sh

EXTRA_DIST += install_pynormaliz.sh

EXTRA_DIST += COPYING INSTALL

# TODO: directory tools? playground? linear_algebra_test?
# TODO: jNormaliz? java libs in lib/? documentation (papers)?

PYNORMALIZ_COMMIT_SHA=f0b1ee424050f25088b77aa4566eb3e33f59b721

dist-hook:
	git clone https://github.com/Normaliz/PyNormaliz $(distdir)/PyNormaliz
	cd $(distdir)/PyNormaliz && git checkout $(PYNORMALIZ_COMMIT_SHA)
	rm -rf $(distdir)/PyNormaliz/.git

