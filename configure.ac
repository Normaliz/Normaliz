dnl Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)

AC_INIT(Normaliz, 3.7.4)

AC_SUBST(LIBNORMALIZ_VERSION_MAJOR, 3)
AC_SUBST(LIBNORMALIZ_VERSION_MINOR, 7)
AC_SUBST(LIBNORMALIZ_VERSION_PATCH, 4)
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

# We can use SCIP (Solving Constraint Integer Programs) for cone decomposition
AC_ARG_WITH([scip], AS_HELP_STRING([--with-scip=FLAGS], [Use SCIP for cone decomposition; the FLAGS specify the linker flags necessary, e.g., -lscip. Typically, a subset of -lscip -lz -lreadline -lscipopt -lnlpi.cppad -llpispx -lsoplex or similar is going to make this work.]))
AS_IF([test "x$with_scip" = "x"], [with_scip=no], [
    AS_IF([test "x$with_scip" = "xyes"], [
          # if no parameters have been given, we try to figure out how to link in SCIP
          AC_SEARCH_LIBS([SCIPversion], [scip scipopt], [], [])
    ], [
          # otherwise, we add the given parameters to LIBS, we'll check below that they actually work.
          LIBS="$LIBS $with_scip"
    ])
])
AS_IF([test "x$with_scip" != "xno"], [
      AC_CHECK_HEADERS([scip/scip.h], [], AC_MSG_ERROR([SCIP headers not found. Did you forget to specify CPPFLAGS=-Idirectory/of/csip/headers?]))
      AC_SEARCH_LIBS([SCIPversion], [], [], AC_MSG_ERROR([SCIP libraries not found. Specify --with-scip='-lscip -lanother -llibrary' if -lscip is not enough to link in SCIP. Did you forget to specify LDFLAGS=-Ldirectory/of/scip/library?]))
      AC_DEFINE(NMZ_SCIP)
      DEFINE_NMZSCIP="#define NMZ_SCIP"
      AC_SUBST([DEFINE_NMZSCIP])
])
AM_CONDITIONAL(NMZ_SCIP, test "x$with_scip" != "xno")

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

AC_CONFIG_FILES([Makefile source/Makefile source/libnormaliz/version.h example/Makefile test/Makefile source/libnormaliz/nmz_config.h])
AC_OUTPUT
