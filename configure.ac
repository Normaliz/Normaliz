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
	example/Makefile
	test/Makefile])
AC_CONFIG_MACRO_DIR([m4])
AC_PROG_CXX
AX_CXX_COMPILE_STDCXX_11(, optional)
AS_IF([test x$HAVE_CXX11 != x1],
  dnl Testing for C++ 11 above is too strict, as
  dnl it rejects the g++ 4.6.3 that is available on Travis.
  [AS_WARN([Did not find a compiler that supports C++ 11.  We keep going with C++0x])
   AX_CXX_COMPILE_STDCXX_0X
])
AC_PROG_LIBTOOL
AC_LANG(C++)

## Test for OpenMP parallelization.
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

## Test whether to build with SCIP.
AC_ARG_WITH([scipoptsuite-src],
   AS_HELP_STRING([--with-scipoptsuite-src=DIR],
      [provide location of scipoptsuite source directory]),
   [
    for a in "$with_scipoptsuite_src"/scip-*.*.*@<:@0-9@:>@ ; do
        AS_IF([test -d $a ],
	   [AS_IF([test x$scipdir != x ],
	       [AC_MSG_ERROR([Directory $with_scipoptsuite_src contains more than one scip directory])])
	    scipdir=$a;
	   ])
    done
    AS_IF([test x$scipdir = x ],
        [AC_MSG_ERROR([Directory $with_scipoptsuite_src does not contain a scip-VERSION directory. Check if it is the right directory and you compiled the scipoptsuite.])])
    AC_MSG_RESULT([Found scipoptsuite in $scipdir])
    scipincludedir="$scipdir"/src
    sciplibdir="$scipdir"/lib
    LDFLAGS="-L$sciplibdir $LDFLAGS"
    CPPFLAGS="-I $scipincludedir $CPPFLAGS"
   ])
AC_ARG_ENABLE([scip],
    [AS_HELP_STRING([--enable-scip@<:@=ARG@:>@],
      [enable cone decomposition with SCIP @<:@default=check@:>@])],
    [:],
    [enable_scip=check])
AS_IF([test "x$enable_scip" != xno],
    [AC_MSG_CHECKING([whether SCIP headers and library are available])
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
    ])
AS_IF([test x$have_scip = xyes],
  [enable_scip=yes],
  [AS_IF([test "x$enable_scip" = xyes],
    [AC_MSG_ERROR([SCIP library is not available])],
    [AS_IF([test "x$enable_scip" != xno],
        [AC_MSG_NOTICE([Cone decomposition with SCIP integer programming solver will not be available because the SCIP library was not found])])
     enable_scip=no])])
AC_SUBST(SCIP_CXXFLAGS)
AC_SUBST(SCIP_LDFLAGS)
AS_IF([test x$enable_scip = xyes],
  [AC_DEFINE(NMZ_SCIP)])
AM_CONDITIONAL(NMZ_SCIP, test x$enable_scip = xyes)

### Test whether to build nmzIntegrate (depends on CoCoALib).
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
     COCOA_LDFLAGS="-lcocoa -lgmp"
     LDFLAGS_SAVED="$LDFLAGS"
     LDFLAGS="$LDFLAGS $COCOA_LDFLAGS"
     export LDFLAGS
     AC_LINK_IFELSE(
       [AC_LANG_PROGRAM([[#include "CoCoA/library.H"
			  using namespace CoCoA;
			 ]],
		    [[   GlobalManager CoCoAFoundations;
		     ]])],
       [have_cocoa=yes],
       [have_cocoa=no])
     AC_MSG_RESULT($have_cocoa)
     LDFLAGS="$LDFLAGS_SAVED"
    ])
AS_IF([test x$have_cocoa = xyes],
  [enable_nmzintegrate=yes],
  [AS_IF([test "x$enable_nmzintegrate" = xyes],
    [AC_MSG_ERROR([CoCoALib is not available and is required for building nmzIntegrate])],
    [AS_IF([test "x$enable_nmzintegrate" != xno],
        [AC_MSG_NOTICE([nmzIntegrate will not be built because CoCoALib was not found])])
     enable_nmzintegrate=no])])
AC_SUBST(COCOA_LDFLAGS)
AC_SUBST(COCOA_CPPFLAGS)
AM_CONDITIONAL(NMZ_INTEGRATE, test x$enable_nmzintegrate = xyes)

AC_OUTPUT
