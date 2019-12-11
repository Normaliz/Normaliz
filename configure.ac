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
   AS_HELP_STRING([--with-gmp@<:@=DIR@:>@],
                  [Use the GMP library installed in installation prefix DIR.]),
   [], [with_gmp=yes])

GMP_CPPFLAGS=
GMP_LDFLAGS=
GMP_LIBS="-lgmpxx -lgmp"
AS_CASE([$with_gmp],
  [no],[
    AC_MSG_ERROR([Building without GMP is not supported])
  ],
  [yes],[],
  [*],[
    GMP_CPPFLAGS="-I${with_gmp}/include"
    GMP_LDFLAGS="-L${with_gmp}/lib"
  ]
)

dnl try to link against GMP
dnl AX_CHECK_LIBRARY([GMP], [gmp.h], [gmp], [__gmpz_init],
dnl     [], AC_MSG_ERROR([Unable to find GMP]))
AX_CHECK_LIBRARY([GMPXX], [gmpxx.h], [gmpxx], [__gmpz_init],
    [], [AC_MSG_ERROR([GMP C++ library not found. Make sure it was compiled with --enable-cxx])], [-lgmp])

CPPFLAGS="$CPPFLAGS $GMP_CPPFLAGS"
LDFLAGS="$LDFLAGS $GMP_LDFLAGS"

AC_MSG_CHECKING([whether a simgple GMP C++ program compiles and links])
LIBS_SAVED="$LIBS"
LIBS="$LIBS $GMP_LIBS"
AC_LINK_IFELSE(
   [AC_LANG_PROGRAM([[#include <gmpxx.h>
                     ]],
                    [[mpz_class a;
                     ]])],
   [AC_MSG_RESULT([yes])],
   [AC_MSG_ERROR([GMP C++ library not found. Make sure it was compiled with --enable-cxx])])
LIBS="$LIBS_SAVED"

AC_SUBST([GMP_LIBS])


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
   AS_HELP_STRING([--with-cocoalib@<:@=DIR@:>@],
      [provide location of CoCoALib installation prefix or source directory]),
   [],[with_cocoalib=check])

COCOA_CPPFLAGS=
COCOA_LDFLAGS=
COCOA_LIBS=
AS_CASE([$with_cocoalib],
  [no],[],
  [yes],[],
  [check],[],
  [*],[
    COCOA_CPPFLAGS="-I$with_cocoalib/include"
    COCOA_LDFLAGS="-L$with_cocoalib/lib"
  ]
)

AS_IF([test "x$with_cocoalib" != xno],
    [AC_MSG_CHECKING([whether CoCoALib headers and library are available])
     COCOA_LIBS="-lcocoa"
     #save_CPPFLAGS="$CPPFLAGS"
     #save_LDFLAGS="$LDFLAGS"
     save_LIBS="$LIBS"
     CPPFLAGS="$CPPFLAGS $COCOA_CPPFLAGS"
     LDFLAGS="$LDFLAGS $COCOA_LDFLAGS"
     LIBS="$LIBS $COCOA_LIBS $GMP_LIBS"
     AC_LINK_IFELSE(
       [AC_LANG_PROGRAM([[#include "CoCoA/library.H"
			  using namespace CoCoA;
			 ]],
		    [[   GlobalManager CoCoAFoundations;
		     ]])],
       [have_cocoa=yes],
       [have_cocoa=no])
     AC_MSG_RESULT($have_cocoa)
     #CPPFLAGS="$save_CPPFLAGS"
     #LDFLAGS="$save_LDFLAGS"
     LIBS="$save_LIBS"
     AS_IF([test x$have_cocoa != xyes],
       [COCOA_LIBS=
        AS_IF([test "x$with_cocoalib" != xcheck],
         [AC_MSG_ERROR([CoCoALib is not available but was requested])])])
    ])

#AC_SUBST(COCOA_CPPFLAGS)
#AC_SUBST(COCOA_LDFLAGS)
AC_SUBST(COCOA_LIBS)
AM_CONDITIONAL(NMZ_COCOA, test x$have_cocoa = xyes)
AS_IF([test x$have_cocoa = xyes],
  [AC_DEFINE(NMZ_COCOA)
   DEFINE_NMZCOCOA="#define NMZ_COCOA"
   AC_SUBST([DEFINE_NMZCOCOA])])


dnl
dnl Test whether to use flint
dnl
AC_ARG_WITH([flint],
   AS_HELP_STRING([--with-flint@<:@=DIR@:>@],
      [provide location of flint installation prefix or source directory]),
   [],[with_flint=check])

FLINT_CPPFLAGS=
FLINT_LDFLAGS=
FLINT_LIBS=
AS_CASE([$with_flint],
  [no],[],
  [yes],[],
  [check],[],
  [*],[
    FLINT_CPPFLAGS="-I$with_flint/include"
    FLINT_LDFLAGS="-L$with_flint/lib -Wl,-rpath,$with_flint/lib"
  ]
)

# TODO/FIXME: mpfr???
AS_IF([test "x$with_flint" != xno],
    [AC_MSG_CHECKING([whether flint headers and library are available])
     FLINT_LIBS="-lflint -lmpfr"
     AX_CHECK_LIBRARY([FLINT], [flint/flint.h], [flint], [fmpz_poly_set_coeff_mpz],[],[],[$GMP_LIBS])
     AS_IF([test x$ax_cv_have_FLINT = xyes],
       [
         CPPFLAGS="$CPPFLAGS $FLINT_CPPFLAGS"
         LDFLAGS="$LDFLAGS $FLINT_LDFLAGS"
       ],
       [FLINT_LIBS=
        AS_IF([test "x$with_flint" != xcheck],
         [AC_MSG_ERROR([flint is not available but was requested])])])
    ])

AC_SUBST(FLINT_LIBS)
AM_CONDITIONAL(NMZ_FLINT, test x$ax_cv_have_FLINT = xyes)
AS_IF([test x$ax_cv_have_FLINT = xyes],
  [AC_DEFINE(NMZ_FLINT)
   DEFINE_NMZFLINT="#define NMZ_FLINT"
   AC_SUBST([DEFINE_NMZFLINT])])


dnl
dnl Test whether to build Normaliz with nauty
dnl
AC_ARG_WITH([nauty],
   AS_HELP_STRING([--with-nauty@<:@=DIR@:>@],
      [provide location of nauty installation prefix or source directory]),
   [],[with_nauty=check])

NAUTY_CPPFLAGS=
NAUTY_LDFLAGS=
NAUTY_LIBS=
AS_CASE([$with_nauty],
  [no],[],
  [yes],[],
  [check],[],
  [*],[
    NAUTY_CPPFLAGS="-I$with_nauty/include"
    NAUTY_LDFLAGS="-L$with_nauty/lib -Wl,-rpath,$with_nauty/lib"
  ]
)

AS_IF([test "x$with_nauty" != xno],
    [AC_MSG_CHECKING([whether nauty headers and library are available])
     NAUTY_LIBS="-lnauty"
     AX_CHECK_LIBRARY([NAUTY], [nauty/nauty.h], [nauty], [densenauty],[],[],[$GMP_LIBS])
     AS_IF([test x$ax_cv_have_NAUTY = xyes],
       [
         CPPFLAGS="$CPPFLAGS $NAUTY_CPPFLAGS"
         LDFLAGS="$LDFLAGS $NAUTY_LDFLAGS"
       ],
       [NAUTY_LIBS=
        AS_IF([test "x$with_nauty" != xcheck],
         [AC_MSG_ERROR([nauty is not available but was requested])])])
    ])

AC_SUBST(NAUTY_LIBS)
AM_CONDITIONAL(NMZ_NAUTY, test x$ax_cv_have_NAUTY = xyes)
AS_IF([test x$ax_cv_have_NAUTY = xyes],
  [AC_DEFINE(NMZ_NAUTY)
   DEFINE_NMZNAUTY="#define NMZ_NAUTY"
   AC_SUBST([DEFINE_NMZNAUTY])])


dnl
dnl Test whether to build Normaliz with e-antic
dnl
AC_ARG_WITH([e-antic],
   AS_HELP_STRING([--with-e-antic@<:@=DIR@:>@],
      [provide location of e-antic installation prefix]),
   [],[with_e_antic=check])

E_ANTIC_CPPFLAGS=
E_ANTIC_LDFLAGS=
E_ANTIC_LIBS=
AS_CASE([$with_e_antic],
  [no],[],
  [yes],[],
  [check],[],
  [*],[
    NAUTY_CPPFLAGS="-I$with_e_antic/include"
    NAUTY_LDFLAGS="-L$with_e_antic/lib -Wl,-rpath,$with_e_antic/lib"
  ]
)

AS_IF([test "x$with_e_antic" != xno],
    [AC_MSG_CHECKING([whether e-antic headers and library are available])
     E_ANTIC_LIBS="-leanticxx -leantic -larb -lflint"
     LIBS_SAVED="$LIBS"
     LIBS="$LIBS $E_ANTIC_LIBS"
     AC_LINK_IFELSE(
       [AC_LANG_PROGRAM([[#include <e-antic/renfxx.h>
                         ]],
                        [[renf_elem_class f;
                         ]])],
       [have_e_antic=yes],
       [have_e_antic=no
        E_ANTIC_LIBS=
        AS_IF([test "x$with_e_antic" != xcheck],
         [AC_MSG_ERROR([e-antic is not available but was requested])])])
     AC_MSG_RESULT($have_e_antic)
     LIBS="$LIBS_SAVED"
    ])
AC_SUBST(E_ANTIC_LIBS)
AM_CONDITIONAL(ENFNORMALIZ, test x$have_e_antic = xyes)
AS_IF([test x$have_e_antic = xyes],
    [AC_DEFINE(ENFNORMALIZ)
    DEFINE_ENFNORMALIZ="#define ENFNORMALIZ"
    AC_SUBST([DEFINE_ENFNORMALIZ])])


dnl
dnl
dnl
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
