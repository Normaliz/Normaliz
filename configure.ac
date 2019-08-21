#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(Normaliz, 3.7.4)
AC_SUBST(LIBNORMALIZ_VERSION_MAJOR, 3)
AC_SUBST(LIBNORMALIZ_VERSION_MINOR, 7)
AC_SUBST(LIBNORMALIZ_VERSION_PATCH, 4)
AC_SUBST(LIBNORMALIZ_VERSION_STRING, "$PACKAGE_VERSION")
AM_INIT_AUTOMAKE(foreign)
AC_CONFIG_FILES([Makefile
	source/Makefile
	source/libnormaliz/version.h
	example/Makefile
	test/Makefile
	source/libnormaliz/nmz_config.h])
AC_CONFIG_MACRO_DIR([m4])
AC_PROG_CXX
AX_CXX_COMPILE_STDCXX(14, , optional)
AS_IF([test x$HAVE_CXX14 = x0],
  [ AX_CXX_COMPILE_STDCXX(0x, , mandatory) ])

AC_PROG_LIBTOOL
AC_LANG(C++)

## Test for GMP
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
AM_CONDITIONAL(NMZ_SCIP, test x$enable_scip = xyes)
AS_IF([test x$enable_scip = xyes],
  [AC_DEFINE(NMZ_SCIP)
   DEFINE_NMZSCIP="#define NMZ_SCIP"
   AC_SUBST([DEFINE_NMZSCIP])])

### Test whether to build Normaliz with CoCoALib.
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

### Test whether to use flint
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
  
### Test whether to build Normaliz with nauty.
AC_ARG_WITH([nautylib],
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

### Test whether to build Normaliz with e-antic
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
     E_ANTIC_LIBS="-leanticxx -leantic -larb -lflint -lmpfr"
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


AC_OUTPUT

