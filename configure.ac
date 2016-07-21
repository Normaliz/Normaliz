#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(normaliz, 3.1.1+autotools-2016-07-21)
AC_SUBST(LIBNORMALIZ_VERSION_MAJOR, 3)
AC_SUBST(LIBNORMALIZ_VERSION_MINOR, 1)
AC_SUBST(LIBNORMALIZ_VERSION_PATCH, 1)
AC_SUBST(LIBNORMALIZ_VERSION_STRING, "$PACKAGE_VERSION")
AM_INIT_AUTOMAKE(foreign)
AC_CONFIG_FILES([Makefile
	source/Makefile
	source/libnormaliz/version.h
	test/Makefile])
AC_CONFIG_MACRO_DIRS([m4])
AC_PROG_CXX
AX_CXX_COMPILE_STDCXX_11(noext, mandatory)
AC_PROG_LIBTOOL
AC_LANG(C++)
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

dnl AC_ARG_WITH([scipoptsuite],
dnl    AS_HELP_STRING([--with-scipoptsuite=DIR],
dnl       [provide location of scipoptsuite source directory]),
dnl    [
dnl      ###  LDFLAGS und CPPFLAGS setzen......
dnl    ]
dnl # with-scipoptsuite=DIR --- for uninstalled SCIP
dnl # with-scip=DIR -- for installed SCIP

dnl #AC_ARG_WITH....

AC_MSG_CHECKING([whether SCIP headers and library are available])
SCIP_LDFLAGS="-lscip -lgmp -lz -lreadline"
LDFLAGS_SAVED="$LDFLAGS"
LDFLAGS="$LDFLAGS $SCIP_LDFLAGS"
export LDFLAGS
AC_LINK_IFELSE(
  [AC_LANG_PROGRAM([[#include <scip/scip.h>
  			    ]],
               [[ SCIP* scip = NULL;
	          SCIPcreate(& scip); ]])],
  [have_scip=yes],
  [have_scip=no])
AC_MSG_RESULT($have_scip)
LDFLAGS="$LDFLAGS_SAVED"

AC_SUBST(SCIP_CXXFLAGS)
AC_SUBST(SCIP_LDFLAGS)
AS_IF([test x$have_scip = xyes],
  [AC_DEFINE(NMZ_SCIP)])
AM_CONDITIONAL(NMZ_SCIP, test x$have_scip = xyes)

AC_OUTPUT
