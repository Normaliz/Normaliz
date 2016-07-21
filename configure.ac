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

AC_OUTPUT
