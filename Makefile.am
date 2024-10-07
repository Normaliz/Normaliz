ACLOCAL_AMFLAGS = -I m4

SUBDIRS = source test example

EXTRA_DIST = bootstrap.sh CHANGELOG \
             .github/workflows/ci.yml \
             .github/workflows/.github-ci-build_configure.sh \
             .github/workflows/.github-ci-build_check.sh \
             .github/workflows/.github-ci-build_make.sh

EXTRA_DIST += Macaulay2/INSTALL Macaulay2/Normaliz.m2
EXTRA_DIST += Singular/README Singular/normaliz.lib Singular/sagbiNormaliz.lib \
              Singular/normaliz.pdf
EXTRA_DIST += doc/Normaliz.pdf
EXTRA_DIST += doc/Normaliz.tex
EXTRA_DIST += doc/Intro.tex
EXTRA_DIST += doc/Discrete.tex
EXTRA_DIST += doc/Monoid.tex
EXTRA_DIST += doc/Input.tex
EXTRA_DIST += doc/Options.tex
EXTRA_DIST += doc/Advanced.tex
EXTRA_DIST += doc/Algebraic.tex
EXTRA_DIST += doc/OptFiles.tex
EXTRA_DIST += doc/Technical.tex
EXTRA_DIST += doc/Appendix.tex
EXTRA_DIST += doc/NmzShortRef.pdf
EXTRA_DIST += doc/NmzShortRef.tex
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
EXTRA_DIST += install_scripts_opt/install_nmz_hash-library.sh
EXTRA_DIST += install_scripts_opt/install_nmz_nauty.sh
EXTRA_DIST += install_scripts_opt/install_nmz_e-antic.sh
EXTRA_DIST += install_scripts_opt/install_eantic_with_prerequisites.sh

EXTRA_DIST += install_scripts_opt/cocoa_patches/eterms.h
EXTRA_DIST += install_scripts_opt/cocoa_patches/SignalWatcher.C
EXTRA_DIST += install_scripts_opt/cocoa_patches/SignalWatcher.H

EXTRA_DIST += install_scripts_opt/e-antic_pataches/nf_elem.h
EXTRA_DIST += install_scripts_opt/e-antic_pataches/nf_elem_add_fmpq.c

EXTRA_DIST += install_scripts_opt/header_files_for_Makefile.classic/nmz_config.h

EXTRA_DIST += source/chunk/run_single.sh

EXTRA_DIST += install_pynormaliz.sh

EXTRA_DIST += COPYING INSTALL

# TODO: directory tools? playground? linear_algebra_test?
# TODO: jNormaliz? java libs in lib/? documentation (papers)?

PYNORMALIZ_COMMIT_SHA=8365df9043492d9f3f83823db6af2b89920706c1

dist-hook:
	git clone https://github.com/Normaliz/PyNormaliz $(distdir)/PyNormaliz
	cd $(distdir)/PyNormaliz && git checkout $(PYNORMALIZ_COMMIT_SHA)
	rm -rf $(distdir)/PyNormaliz/.git


