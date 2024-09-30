dnl
dnl Process this file with autoconf to produce a configure script.
dnl


dnl
dnl Init autoconf, set package version
dnl
m4_define([normaliz_major_version], [3])
m4_define([normaliz_minor_version], [10])
m4_define([normaliz_patch_version], [4])
m4_define([normaliz_version],
      [normaliz_major_version.normaliz_minor_version.normaliz_patch_version])

AC_PREREQ([2.68])
AC_INIT([Normaliz],
        [normaliz_version],
        [https://github.com/Normaliz/Normaliz/issues],
        [normaliz],
        [https://www.normaliz.uni-osnabrueck.de])


AC_SUBST(LIBNORMALIZ_VERSION_MAJOR, normaliz_major_version)
AC_SUBST(LIBNORMALIZ_VERSION_MINOR, normaliz_minor_version)
AC_SUBST(LIBNORMALIZ_VERSION_PATCH, normaliz_patch_version)
AC_SUBST(LIBNORMALIZ_VERSION_STRING, "$PACKAGE_VERSION")

AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_AUX_DIR([cnf])


dnl
dnl Init automake and libtool
dnl
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
CPPFLAGS="$CPPFLAGS $GMP_CPPFLAGS"
LDFLAGS="$LDFLAGS $GMP_LDFLAGS"
dnl AX_CHECK_LIBRARY([GMP], [gmp.h], [gmp], [__gmpz_init],
dnl     [], AC_MSG_ERROR([Unable to find GMP]))
AX_CHECK_LIBRARY([GMPXX], [gmpxx.h], [gmpxx], [__gmpz_init],
    [], [AC_MSG_ERROR([GMP C++ library not found. Make sure it was compiled with --enable-cxx])], [-lgmp])

AC_MSG_CHECKING([whether a simple GMP C++ program compiles and links])
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
    FLINT_LDFLAGS="-L$with_flint/lib"
  ]
)

# TODO/FIXME: mpfr???
AS_IF([test "x$with_flint" != xno],
    [AC_MSG_CHECKING([whether flint headers and library are available])
     FLINT_LIBS="-lflint -lmpfr"
     AX_CHECK_LIBRARY([FLINT], [flint/flint.h], [flint], [fmpz_poly_set_coeff_fmpz],[],[],[$GMP_LIBS])
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
dnl Test whether to build Normaliz with hash-library
dnl
AC_ARG_WITH([hashlibrary],
   AS_HELP_STRING([--with-hashlibrary@<:@=DIR@:>@],
      [provide location of hash-library installation prefix or source directory]),
   [],[with_hashlibrary=check])

HASHLIBRARY_CPPFLAGS=
HASHLIBRARY_LDFLAGS=
HASHLIBRARY_LIBS=
AS_CASE([$with_hashlibrary],
  [no],[],
  [yes],[],
  [check],[],
  [*],[
    HASHLIBRARY_CPPFLAGS="-I$with_hashlibrary/include"
    HASHLIBRARY_LDFLAGS="-L$with_hashlibrary/lib"
  ]
)

AS_IF([test "x$with_hashlibrary" != xno],
    [AC_MSG_CHECKING([whether hash-library headers and library are available])
     HASHLIBRARY_LIBS="-lsha256"
     #save_CPPFLAGS="$CPPFLAGS"
     #save_LDFLAGS="$LDFLAGS"
     save_LIBS="$LIBS"
     LIBS="$LIBS $HASHLIBRARY_LIBS"
     AC_LINK_IFELSE(
       [AC_LANG_PROGRAM([[#include "hash-library/sha256.h"
                        ]],
                        [[SHA256 dummy;
                        ]])],
       [[have_hashlibrary=yes]
        [CPPFLAGS="$CPPFLAGS $HASHLIBRARY_CPPFLAGS"]
        [LDFLAGS="$LDFLAGS $HASHLIBRARY_LDFLAGS"]],
       [[have_hashlibrary=no]
        [LIBS="$save_LIBS"]])
     AC_MSG_RESULT($have_hashlibrary)
     # LIBS="$save_LIBS"
     AS_IF([test x$have_hashlibrary != xyes],
       [HASHLIBRARY_LIBS=
        AS_IF([test "x$with_hashlibrary" != xcheck],
         [AC_MSG_ERROR([hash-library is not available but was requested])])])
    ])

AC_SUBST(HASHLIBRARY_LIBS)
AM_CONDITIONAL(NMZ_HASHLIBRARY, test x$have_hashlibrary = xyes)
AS_IF([test x$have_hashlibrary = xyes],
  [AC_DEFINE(NMZ_HASHLIBRARY)
   DEFINE_NMZHASHLIBRARY="#define NMZ_HASHLIBRARY"
   AC_SUBST([DEFINE_NMZHASHLIBRARY])])


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
    NAUTY_LDFLAGS="-L$with_nauty/lib"
  ]
)

dnl
dnl Check whether nauty header is "nauty/nauty.h".
dnl If not, check for "nauty.h" below.
dnl
AS_IF([test "x$with_nauty" != xno],
    [AC_MSG_CHECKING([whether nauty headers and library are available (nauty/nauty.h)])
     NAUTY_LIBS="-lnauty"
     AX_CHECK_LIBRARY([NAUTYNAUTY], [nauty/nauty.h], [nauty], [densenauty],[],[],[$GMP_LIBS])
     AS_IF([test x$ax_cv_have_NAUTYNAUTY = xyes],
       [
         CPPFLAGS="$CPPFLAGS $NAUTY_CPPFLAGS"
         LDFLAGS="$LDFLAGS $NAUTY_LDFLAGS"
         AC_MSG_CHECKING([whether nauty supports TLS])
         AC_COMPILE_IFELSE(
           [AC_LANG_PROGRAM([[#include "nauty/nauty.h"
                            ]],
                            [[#if(!(defined HAVE_TLS) || 0 == HAVE_TLS)
                                #error "nauty without TLS according to nauty/nauty.h"
                              #endif
                            ]])],
           [nauty_have_tls=yes],
           [nauty_have_tls=no])
         AC_MSG_RESULT($nauty_have_tls)
       ],
       [NAUTY_LIBS=
       ])
    ])

dnl
dnl If "nauty/nauty.h" not found above, check for "nauty.h".
dnl
AS_IF([test "x$with_nauty" != xno && test "x$ax_cv_have_NAUTYNAUTY" != xyes],
    [AC_MSG_CHECKING([whether nauty headers and library are available (nauty.h)])
     NAUTY_LIBS="-lnauty"
     AX_CHECK_LIBRARY([NAUTY], [nauty.h], [nauty], [densenauty],[],[],[$GMP_LIBS])
     AS_IF([test x$ax_cv_have_NAUTY = xyes],
       [
         CPPFLAGS="$CPPFLAGS $NAUTY_CPPFLAGS"
         LDFLAGS="$LDFLAGS $NAUTY_LDFLAGS"
         AC_MSG_CHECKING([whether nauty supports TLS])
         AC_COMPILE_IFELSE(
           [AC_LANG_PROGRAM([[#include "nauty.h"
                            ]],
                            [[#if(!(defined HAVE_TLS) || 0 == HAVE_TLS)
                                #error "nauty without TLS according to nauty.h"
                              #endif
                            ]])],
           [nauty_have_tls=yes],
           [nauty_have_tls=no])
         AC_MSG_RESULT($nauty_have_tls)
       ],
       [NAUTY_LIBS=
       ])
    ])

AS_IF([test "x$with_nauty" != xno && test "x$with_nauty" != xcheck && test "x$ax_cv_have_NAUTYNAUTY" != xyes && test "x$ax_cv_have_NAUTY" != xyes],
      [AC_MSG_ERROR([nauty is not available but was requested])])

AC_SUBST(NAUTY_LIBS)
AM_CONDITIONAL(NMZ_NAUTYNAUTY, [test x$ax_cv_have_NAUTYNAUTY = xyes])
AM_CONDITIONAL(NMZ_NAUTY, [test x$ax_cv_have_NAUTYNAUTY = xyes || test x$ax_cv_have_NAUTY = xyes])
AM_CONDITIONAL(NMZ_NAUTY_TLS, [test x$nauty_have_tls = xyes])

AS_IF([test x$ax_cv_have_NAUTYNAUTY = xyes],
  [
   AC_DEFINE(NMZ_NAUTYNAUTY)
   DEFINE_NMZNAUTYNAUTY="#define NMZ_NAUTYNAUTY"
   AC_SUBST([DEFINE_NMZNAUTYNAUTY])
   AC_DEFINE(NMZ_NAUTY)
   DEFINE_NMZNAUTY="#define NMZ_NAUTY"
   AC_SUBST([DEFINE_NMZNAUTY])
   AS_IF([test x$nauty_have_tls = xyes],
     [AC_DEFINE(NMZ_NAUTY_TLS)
      DEFINE_NMZNAUTYTLS="#define NMZ_NAUTY_TLS"
      AC_SUBST([DEFINE_NMZNAUTYTLS])
     ])
  ])

AS_IF([test x$ax_cv_have_NAUTY = xyes],
  [
   AC_DEFINE(NMZ_NAUTY)
   DEFINE_NMZNAUTY="#define NMZ_NAUTY"
   AC_SUBST([DEFINE_NMZNAUTY])
   AS_IF([test x$nauty_have_tls = xyes],
     [AC_DEFINE(NMZ_NAUTY_TLS)
      DEFINE_NMZNAUTYTLS="#define NMZ_NAUTY_TLS"
      AC_SUBST([DEFINE_NMZNAUTYTLS])
     ])
  ])


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
    E_ANTIC_CPPFLAGS="-I$with_e_antic/include"
    E_ANTIC_LDFLAGS="-L$with_e_antic/lib"
  ]
)

AS_IF([test "x$with_e_antic" != xno],
    [
     E_ANTIC_LIBS="-leanticxx -leantic"
     AC_MSG_CHECKING([whether e-antic headers and library are available])
     CPPFLAGS_SAVED="$CPPFLAGS"
     LDFLAGS_SAVED="$LDFLAGS"
     LIBS_SAVED="$LIBS"
     CPPFLAGS="$CPPFLAGS $E_ANTIC_CPPFLAGS"
     LDFLAGS="$LDFLAGS $E_ANTIC_LDFLAGS"
     LIBS="$LIBS $E_ANTIC_LIBS"
     AC_LINK_IFELSE(
       [AC_LANG_PROGRAM([[#include <e-antic/renfxx.h>
                         ]],
                        [[ eantic::renf_elem_class f;
                         ]])],
       [have_e_antic=yes],
       [have_e_antic=no
        E_ANTIC_LIBS=
        AS_IF([test "x$with_e_antic" != xcheck],
         [AC_MSG_ERROR([e-antic is not available in a supported version but was requested])])])
     AC_MSG_RESULT($have_e_antic)
     CPPFLAGS="$CPPFLAGS_SAVED"
     LDFLAGS="$LDFLAGS_SAVED"
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
