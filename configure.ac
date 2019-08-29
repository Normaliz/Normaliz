dnl Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)

AC_INIT(Normaliz, 3.8.0)

AC_SUBST(LIBNORMALIZ_VERSION_MAJOR, 3)
AC_SUBST(LIBNORMALIZ_VERSION_MINOR, 8)
AC_SUBST(LIBNORMALIZ_VERSION_PATCH, 0)
AC_SUBST(LIBNORMALIZ_VERSION_STRING, "$PACKAGE_VERSION")

AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_AUX_DIR([.])
AM_INIT_AUTOMAKE([1.14.1])

AC_CANONICAL_HOST

dnl Find C++ Compiler
AC_PROG_CXX

dnl Set INSTALL or use install-sh.
AC_PROG_INSTALL

dnl We build our library with libtool.
LT_INIT

dnl We want to use some C++14 features
AX_CXX_COMPILE_STDCXX(14, , optional)

AS_IF([test x$HAVE_CXX14 = x0],
  dnl We require C++0x features
  [ AX_CXX_COMPILE_STDCXX(0x, , mandatory) ])

dnl Check for required/optional libraries
AC_LANG(C++)

AC_CHECK_HEADERS([boost/dynamic_bitset.hpp], , AC_MSG_ERROR([boost headers not found]))

# We use GMP's C++ wrapper but we also interface with GMP's C interface directly.
AC_CHECK_LIB([gmp], [__gmpz_init], [have_gmp=yes], AC_MSG_ERROR([GMP library not found]))

# GMPXX does not contain anything that we can check for with AX_CXX_CHECK_LIB
# so we just check for something from the standard library, i.e., that -lgmpxx
# goes through
AX_CXX_CHECK_LIB([gmpxx], [std::runtime_error::what () const], [have_gmpxx=yes], AC_MSG_ERROR([GMP library for C++ not found]), [-lgmp])

# We support OpenMP for parallelization if present
AC_ARG_ENABLE([openmp], AS_HELP_STRING([--disable-openmp], [Do not build with support or OpenMP]))
AS_IF([test "x$enable_openmp" != "xno"],
      [
       enable_openmp=yes
       AX_OPENMP(
        [AC_SUBST(OPENMP_CXXFLAGS)],
        AC_MSG_ERROR([We do not know how to enable OpenMP support for your compiler. You can either disable OpenMP with --disable-openmp or set OPENMP_CXXFLAGS to the flags that your compiler understands.]))
      ], [])
AM_CONDITIONAL(ENABLE_OPENMP, [test "x$enable_openmp" = "xyes"])

<<<<<<< HEAD
# We can optionally build with CoCoALib
AC_ARG_WITH([cocoalib], AS_HELP_STRING([--with-cocoalib], [Use CoCoALib for nmzIntegrate.]))
AS_IF([test "x$with_cocoalib" = "x"], [with_cocoalib=no])
AS_IF([test "x$with_cocoalib" != "xno"], [
    AC_CHECK_HEADERS([CoCoA/libray.H], [], AC_MSG_ERROR([CoCoALib headers not found."]))
    AX_CXX_CHECK_LIB([cocoa], [CoCoA::GlobalManager], [], AC_MSG_ERROR([CoCoALib library not found."]))
    AC_DEFINE(NMZ_COCOA)
    DEFINE_NMZCOCOA="#define NMZ_COCOA"
    AC_SUBST([DEFINE_NMZCOCOA])
])
AM_CONDITIONAL(NMZ_INTEGRATE, test "x$with_cocoalib" != "xno")
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

# We can optionally build with FLINT
AC_ARG_WITH([flint], AS_HELP_STRING([--with-flint], [Use FLINT]))
AS_IF([test "x$with_flint" = "x"], [with_flint=no])
AS_IF([test "x$with_flint" != "xno"], [
    AC_CHECK_LIB([flint], [fmpz_init], [], [AC_MSG_ERROR([libflint not found.])], [-lmpfr])
    LIBS="-lmpfr $LIBS"
    AC_DEFINE(NMZ_FLINT)
    DEFINE_NMZFLINT="#define NMZ_FLINT"
    AC_SUBST([DEFINE_NMZFLINT])
])
AM_CONDITIONAL(NMZ_FLINT, test "x$with_flint" != "xno")
  
# We can optionally build with nauty.
AC_ARG_WITH([nauty], AS_HELP_STRING([--with-nauty], [Use nauty for nmzInautomorph]))
AS_IF([test "x$with_nauty" = "x"], [with_nauty=no])
AS_IF([test "x$with_nauty" != "xno"], [
    AC_CHECK_LIB([nauty], [densenauty], [], [AC_MSG_ERROR([nauty library not found.])])
    AC_DEFINE(NMZ_NAUTY)
    DEFINE_NMZNAUTY="#define NMZ_NAUTY"
    AC_SUBST([DEFINE_NMZNAUTY])
])
AM_CONDITIONAL(NMZ_NAUTY, test "x$with_nauty" != "xno")

# We can optionally build with  e-antic
AC_ARG_WITH([e-antic], AS_HELP_STRING([--with-e-antic], [Use e-antic for ENFnormaliz]))
AS_IF([test "x$with_e_antic" = "x"], [with_e_antic=no])
AS_IF([test "x$with_e_antic" != "xno"], [
    AC_CHECK_HEADERS([e-antic/renfxx.h], [], AC_MSG_ERROR([e-antic headers not found]))
    AX_CXX_CHECK_LIB([eanticxx], [eantic::renf_elem_class], [], AC_MSG_ERROR([e-antic library not found]), [-lmpfr -lflint])
    LIBS="-lmpfr -lflint $LIBS"
    AC_DEFINE(ENFNORMALIZ)
    DEFINE_ENFNORMALIZ="#define ENFNORMALIZ"
    AC_SUBST([DEFINE_ENFNORMALIZ])
])
AM_CONDITIONAL(ENFNORMALIZ, test "x$with_e_antic" != "xno")

AS_IF([test x$enable_shared = xyes],
    [AC_DEFINE(DEFINE_NORMALIZ_DLL)
    DEFINE_NORMALIZ_DLL="#define NORMALIZ_USE_DLL"
    AC_SUBST([DEFINE_NORMALIZ_DLL])])

AC_CONFIG_FILES([Makefile source/Makefile source/libnormaliz/version.h example/Makefile test/Makefile source/libnormaliz/nmz_config.h])

AC_OUTPUT
