dnl
dnl Process this file with autoconf to produce a configure script.
dnl


dnl
dnl Init autoconf, set package version
dnl
AC_PREREQ([2.68])
AC_INIT([Normaliz],
        [3.8.3],
        [https://github.com/Normaliz/Normaliz/issues],
        [normaliz],
        [https://www.normaliz.uni-osnabrueck.de])


AC_SUBST(LIBNORMALIZ_VERSION_MAJOR, 3)
AC_SUBST(LIBNORMALIZ_VERSION_MINOR, 8)
AC_SUBST(LIBNORMALIZ_VERSION_PATCH, 3)
AC_SUBST(LIBNORMALIZ_VERSION_STRING, "$PACKAGE_VERSION")

AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_AUX_DIR([cnf])


dnl
dnl Init automake and libtool
AM_INIT_AUTOMAKE([foreign dist-zip subdir-objects std-options])
LT_INIT


dnl
dnl Check for working C++ compiler; ask for C++14, require C++x0
dnl
AC_PROG_CXX
AX_CXX_COMPILE_STDCXX(14, , optional)
AS_IF([test x$HAVE_CXX14 = x0],
  [ AX_CXX_COMPILE_STDCXX(0x, , mandatory) ])

AC_LANG(C++)


dnl
dnl compiler builtins
dnl
AC_DEFUN([CHECK_COMPILER_BUILTIN],
[AC_MSG_CHECKING([for $1])
    AC_LINK_IFELSE(
        [AC_LANG_PROGRAM(
            [[]],
            [$1[($2)];
            ]
        )],
        [AS_VAR_SET([[have_]$1], [yes])],
        [AS_VAR_SET([[have_]$1], [no])]
        )
    AC_MSG_RESULT(AS_VAR_GET([[have_]$1]))
    AS_IF([test yes = AS_VAR_GET([[have_]$1])],
        [AC_DEFINE_UNQUOTED(AS_TR_CPP([HAVE_]$1), 1,
            [Define to 1 if the system has the `]$1[' built-in function])], []
        )])

CHECK_COMPILER_BUILTIN([__builtin_ctzll],[0]);
CHECK_COMPILER_BUILTIN([__builtin_popcountll],[0]);

dnl determine sizeof for some standard types
AC_CHECK_SIZEOF([long long])


dnl
dnl Test for GMP
dnl
AC_ARG_WITH(gmp,
   AS_HELP_STRING([--with-gmp=DIR],
                  [Use the GMP library installed in installation prefix DIR.]),
   [ LDFLAGS="-L$with_gmp/lib $LDFLAGS"
     CPPFLAGS="-I $with_gmp/include $CPPFLAGS"
   ])
AC_MSG_CHECKING([whether GMP C++ library is available])
LIBS_SAVED="$LIBS"
LIBS="$LIBS -lgmpxx -lgmp"
AC_LINK_IFELSE(
   [AC_LANG_PROGRAM([[#include <gmpxx.h>
                     ]],
                    [[mpz_class a;
                     ]])],
   [AC_MSG_RESULT([yes])],
   [AC_MSG_ERROR([GMP C++ library not found. Make sure it was compiled with --enable-cxx])])
LIBS="$LIBS_SAVED"


dnl
dnl Test for OpenMP parallelization
dnl
AC_ARG_ENABLE([openmp],
    [AS_HELP_STRING([--enable-openmp@<:@=ARG@:>@],
      [enable parallelization with OpenMP @<:@default=check@:>@])],
    [:],
    [enable_openmp=check])
AS_IF([test "x$enable_openmp" != xno],
  [AX_OPENMP([have_openmp=yes],
	     [have_openmp=no])])

AS_IF([test x$have_openmp = xyes],
  [enable_openmp=yes],
  [AS_IF([test "x$enable_openmp" = xyes],
    [AC_MSG_ERROR([openmp not available])],
    [AS_IF([test "x$enable_openmp" != xno],
        [AC_MSG_NOTICE([Parallelization with OpenMP will not be available because no suitable compiler was found])])
     enable_openmp=no])])

AC_SUBST(OPENMP_CXXFLAGS)
AM_CONDITIONAL(ENABLE_OPENMP, [test x$enable_openmp = xyes])


dnl
dnl Test whether to build Normaliz with CoCoALib
dnl
AC_ARG_WITH([cocoalib],
   AS_HELP_STRING([--with-cocoalib=DIR],
      [provide location of CoCoALib installation prefix or source directory]),
   [ LDFLAGS="-L$with_cocoalib/lib $LDFLAGS"
     CPPFLAGS="-I $with_cocoalib/include $CPPFLAGS"
   ])
AC_ARG_ENABLE([nmzintegrate],
    [AS_HELP_STRING([--enable-nmzintegrate@<:@=ARG@:>@],
      [enable nmzIntegrate (requires CoCoALib) @<:@default=check@:>@])],
    [:],
    [enable_nmzintegrate=check])
AS_IF([test "x$enable_nmzintegrate" != xno],
    [AC_MSG_CHECKING([whether CoCoALib headers and library are available])
     COCOA_LIBS="-lcocoa -lgmp"
     LIBS_SAVED="$LIBS"
     LIBS="$LIBS $COCOA_LIBS"
     AC_LINK_IFELSE(
       [AC_LANG_PROGRAM([[#include "CoCoA/library.H"
			  using namespace CoCoA;
			 ]],
		    [[   GlobalManager CoCoAFoundations;
		     ]])],
       [have_cocoa=yes],
       [have_cocoa=no])
     AC_MSG_RESULT($have_cocoa)
     LIBS="$LIBS_SAVED"
    ])
AS_IF([test x$have_cocoa = xyes],
  [enable_nmzintegrate=yes],
  [AS_IF([test "x$enable_nmzintegrate" = xyes],
    [AC_MSG_ERROR([CoCoALib is not available and is required for building nmzIntegrate])],
    [AS_IF([test "x$enable_nmzintegrate" != xno],
        [AC_MSG_NOTICE([Normaliz will be built without CoCoALib because it was not found])])
     enable_nmzintegrate=no])])
AC_SUBST(COCOA_LIBS)
AC_SUBST(COCOA_CPPFLAGS)
AM_CONDITIONAL(NMZ_INTEGRATE, test x$enable_nmzintegrate = xyes)
AS_IF([test x$enable_nmzintegrate = xyes],
  [AC_DEFINE(NMZ_COCOA)
   DEFINE_NMZCOCOA="#define NMZ_COCOA"
   AC_SUBST([DEFINE_NMZCOCOA])])


dnl
dnl Test whether to use flint
dnl
AC_ARG_WITH([flint],
   AS_HELP_STRING([--with-flint=DIR],
      [provide location of flint installation prefix or source directory]),
   [ LDFLAGS="-L$with_flint/lib -Wl,-rpath,$with_flint/lib $LDFLAGS"
     CPPFLAGS="-I$with_flint/include $CPPFLAGS"
   ])
AC_ARG_ENABLE([flint],
    [AS_HELP_STRING([--enable-flint@<:@=ARG@:>@],
      [enable flint @<:@default=check@:>@])],
    [:],
    [enable_flint=check])
AS_IF([test "x$enable_flint" != xno],
    [FLINT_LIBS="-lflint -lmpfr"
     LIBS_SAVED="$LIBS"
     LIBS="$LIBS $FLINT_LIBS"
     AC_CHECK_LIB(flint,dummy_fn_to_prevent_ISOC_warning,[have_flint=yes],[have_flint=no])
     LIBS="$LIBS_SAVED"
    ])
AS_IF([test x$have_flint = xyes],
  [enable_flint=yes],
  [AS_IF([test "x$enable_flint" = xyes],
    [AC_MSG_ERROR([flint is not available])],
    [AS_IF([test "x$enable_flint" != xno],
        [AC_MSG_NOTICE([flint will not be used because flint could not be found])])
     enable_flint=no])])
AC_SUBST(FLINT_LIBS)
AC_SUBST(COCOA_CPPFLAGS)
AM_CONDITIONAL(NMZ_FLINT, test "x$enable_flint" = xyes)
AS_IF([test "x$enable_flint" = xyes],
  [AC_DEFINE(NMZ_FLINT)
   DEFINE_NMZFLINT="#define NMZ_FLINT"
   AC_SUBST([DEFINE_NMZFLINT])])


dnl
dnl Test whether to build Normaliz with nauty
dnl
AC_ARG_WITH([nauty],
   AS_HELP_STRING([--with-nauty=DIR],
      [provide location of nauty installation prefix or source directory]),
   [ LDFLAGS="-L$with_nauty/lib $LDFLAGS"
     CPPFLAGS="-I $with_nauty/include $CPPFLAGS"
   ])
AC_ARG_ENABLE([nmzinzautomorph],
    [AS_HELP_STRING([--enable-nmzautomorph@<:@=ARG@:>@],
      [enable nmzInautomorph (requires nauty) @<:@default=check@:>@])],
    [:],
    [enable_nmznauty=check])
AS_IF([test "x$enable_nmznauty" != xno],
    [AC_MSG_CHECKING([whether nauty headers and library are available])
     NAUTY_LIBS="-lnauty"
     LIBS_SAVED="$LIBS"
     LIBS="$LIBS $NAUTY_LIBS"
     AC_CHECK_LIB(nauty,densenauty,[have_nauty=yes],[have_nauty=no])
     AC_MSG_RESULT($have_nauty)
     LIBS="$LIBS_SAVED"
    ])
AS_IF([test x$have_nauty = xyes],
  [enable_nmznauty=yes],
  [AS_IF([test "x$enable_nmznauty" = xyes],
    [AC_MSG_ERROR([nauty is not available and is required for automorphisms])],
    [AS_IF([test "x$enable_nmzintegrate" != xno],
        [AC_MSG_NOTICE([normaliz builzauztomorphisms  built because nauty was not found])])
     enable_nmznauty=no])])
AC_SUBST(NAUTY_LIBS)
AC_SUBST(NAITY_CPPFLAGS)
AM_CONDITIONAL(NMZ_NAUTY, test x$enable_nmznauty = xyes)
AS_IF([test "x$have_nauty" = xyes],
  [AC_DEFINE(NMZ_NAUTY)
   DEFINE_NMZNAUTY="#define NMZ_NAUTY"
   AC_SUBST([DEFINE_NMZNAUTY])])


dnl
dnl Test whether to build Normaliz with e-antic
dnl
AC_ARG_WITH([e-antic],
   AS_HELP_STRING([--with-e-antic=DIR],
      [provide location of e-antic installation prefix]),
   [ LDFLAGS="-L$with_e_antic/lib $LDFLAGS"
     CPPFLAGS="-I $with_e_antic/include $CPPFLAGS"
   ])
AC_ARG_ENABLE([enfnormaliz],
    [AS_HELP_STRING([--enable-enfnormaliz@<:@=ARG@:>@],
      [enable ENFnormaliz (requires e-antic) @<:@default=check@:>@])],
    [:],
    [enable_enfnormaliz=check])
AS_IF([test "x$enable_enfnormaliz" != xno],
    [AC_MSG_CHECKING([whether e-antic headers and library are available])
     E_ANTIC_LIBS="-leanticxx -leantic -larb -lflint"
     LIBS_SAVED="$LIBS"
     LIBS="$LIBS $E_ANTIC_LIBS"
     AC_LINK_IFELSE(
       [AC_LANG_PROGRAM([[#include <e-antic/renfxx.h>
                         ]],
                        [[renf_elem_class f;
                         ]])],
       [have_e_antic=yes; enable_enfnormaliz=yes],
       [have_e_antic=no; enable_enfnormaliz=no])
     AC_MSG_RESULT($have_e_antic)
     LIBS="$LIBS_SAVED"
    ])
AC_SUBST(E_ANTIC_LIBS)
AM_CONDITIONAL(ENFNORMALIZ, test x$enable_enfnormaliz = xyes)
AS_IF([test x$enable_enfnormaliz = xyes],
    [AC_DEFINE(ENFNORMALIZ)
    DEFINE_ENFNORMALIZ="#define ENFNORMALIZ"
    AC_SUBST([DEFINE_ENFNORMALIZ])])

AS_IF([test x$enable_shared = xyes],
    [AC_DEFINE(DEFINE_NORMALIZ_DLL)
    DEFINE_NORMALIZ_DLL="#define NORMALIZ_USE_DLL"
    AC_SUBST([DEFINE_NORMALIZ_DLL])])


dnl
dnl Output everything
dnl
AC_CONFIG_FILES([
    Makefile
	example/Makefile
	source/libnormaliz/nmz_config.h
	source/libnormaliz/version.h
	source/Makefile
	test/Makefile
	])
AC_OUTPUT
