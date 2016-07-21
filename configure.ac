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

AC_OUTPUT
