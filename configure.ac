#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(normaliz, 3.1.1+autotools-2016-07-23)
AC_SUBST(LIBNORMALIZ_VERSION_MAJOR, 3)
AC_SUBST(LIBNORMALIZ_VERSION_MINOR, 1)
AC_SUBST(LIBNORMALIZ_VERSION_PATCH, 2)
AC_SUBST(LIBNORMALIZ_VERSION_STRING, "$PACKAGE_VERSION")
AM_INIT_AUTOMAKE(foreign)
AC_CONFIG_FILES([Makefile
	source/Makefile
	source/libnormaliz/version.h
	example/Makefile
	test/Makefile])
AC_CONFIG_MACRO_DIR([m4])
AC_PROG_CXX
AX_CXX_COMPILE_STDCXX(0x, , mandatory)

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
    AC_MSG_CHECKING([contents of scipoptsuite source directory])
    # The following tests correspond to scipoptsuite 3.2.0, 3.2.1 (hand-written Makefiles)
    for a in "$with_scipoptsuite_src"/scip-*.*.*@<:@0-9@:>@ ; do
        AS_IF([test -d $a ],
	   [AS_IF([test x$scipdir != x ],
	       [AC_MSG_ERROR([Directory $with_scipoptsuite_src contains more than one scip directory])])
	    scipdir=`cd $a && pwd`
	    scipoptsuitedir=`cd $with_scipoptsuite_src && pwd`
	   ])
    done
    AS_IF([test x$scipdir = x ],
        [AC_MSG_ERROR([Directory $with_scipoptsuite_src does not contain a scip-VERSION directory. Check if it is the right directory and you compiled the scipoptsuite.])])
    AC_MSG_RESULT([Found scipoptsuite in $scipdir])
    scipincludedir="$scipdir"/src
    scipoptsuitelibdir="$scipoptsuitedir"/lib  # this is where libscipopt is built
    sciplibdir="$scipdir"/lib		       # this is where libscip is built
    AS_IF([test -f "$scipoptsuitelibdir"/libscipopt.a ],
        [AC_MSG_RESULT([ Found short link $a ]) ],
        [# short link missing (not installed by 3.2.0, exists in 3.2.1)
	 for a in "$scipoptsuitelibdir"/libscipopt-*.*.a ; do
	    AS_IF([test -f $a ],
	       [AC_MSG_RESULT([Found library $a])
	        AS_IF([test x$scipsuffix != x ],
		   [AC_MSG_ERROR([Directory $scipoptsuitelibdir contains more than one libscipopt*.a library, please create a link to libscipopt.a ])])
		scipsuffix=`echo $a | sed 's|.*/libscipopt\(@<:@^/@:>@*\)\.a|\1|'`
	       ])
         done
	 AS_IF([test x$scipsuffix = x],
	     [AC_MSG_ERROR([Directory $scipoptsuitelibdir does not contain a libscipopt*.a library])])
	])
    LDFLAGS="-L$scipoptsuitelibdir $LDFLAGS"
    CPPFLAGS="-I $scipincludedir $CPPFLAGS"
   ])
AC_ARG_ENABLE([scip],
    [AS_HELP_STRING([--enable-scip@<:@=ARG@:>@],
      [enable cone decomposition with SCIP @<:@default=check@:>@])],
    [:],
    [enable_scip=check])
AS_IF([test "x$enable_scip" != xno],
    [AC_MSG_CHECKING([whether SCIP headers and library are available])
     dnl SCIP can be installed in many different ways.  We try several.
     dnl There's libscipopt and there's libscip (which depends on
     dnl liblpispx and libsoplex or another LP solver library; but these
     dnl dependencies may not be recorded).
     LIBS_SAVED="$LIBS"
     for SCIP_LIBS in \
	 "-lscip" \
	 "-lscipopt$scipsuffix -lgmp -lz -lreadline" \
	 "-lscipopt$scipsuffix -lz -lreadline" \
	 "-lscipopt$scipsuffix -lgmp -lz" \
	 "-lscipopt$scipsuffix -lz" \
	 "-lscipopt$scipsuffix -lgmp -lreadline" \
	 "-lscipopt$scipsuffix -lreadline" \
	 "-lscipopt$scipsuffix -lgmp" \
	 "-lscipopt" \
	 "-lscip -lnlpi.cppad -llpispx -lsoplex -lgmp -lz -lreadline" \
	 ; do
	 LIBS="$LIBS_SAVED $SCIP_LIBS"
	 dnl We use AC_RUN_IF_ELSE because on Mac OS X, running
	 dnl the program is the only way to know whether we got
	 dnl all dylibs working.
	 AC_RUN_IFELSE(
	   [AC_LANG_PROGRAM([[#include <scip/scip.h>
				     ]],
			[[ SCIP* scip = NULL;
			   SCIPcreate(& scip);
			   return 0;
			 ]])],
	   [have_scip=yes
	    AC_MSG_RESULT([yes, using $SCIP_LIBS])
	    break
	   ],
	   [have_scip=no])
     done
     LIBS="$LIBS_SAVED"
     if test $have_scip != yes ; then
	AC_MSG_RESULT($have_scip)
     fi
    ])
AS_IF([test x$have_scip = xyes],
  [enable_scip=yes],
  [AS_IF([test "x$enable_scip" = xyes],
    [AC_MSG_ERROR([SCIP library is not available])],
    [AS_IF([test "x$enable_scip" != xno],
        [AC_MSG_NOTICE([Cone decomposition with SCIP integer programming solver will not be available because the SCIP library was not found])])
     enable_scip=no])])
AC_SUBST(SCIP_CXXFLAGS)
AC_SUBST(SCIP_LIBS)
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
        [AC_MSG_NOTICE([nmzIntegrate will not be built because CoCoALib was not found])])
     enable_nmzintegrate=no])])
AC_SUBST(COCOA_LIBS)
AC_SUBST(COCOA_CPPFLAGS)
AM_CONDITIONAL(NMZ_INTEGRATE, test x$enable_nmzintegrate = xyes)

AC_OUTPUT
