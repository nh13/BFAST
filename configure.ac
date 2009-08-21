#  See LICENSE for details
#
#  The version number in AC_INIT is for the bfast package.
#  Developers - the version number should *NEVER* be changed except by, 
#  or with explicit permission from, Nils Homer.

AC_PREREQ(2.59)
AC_INIT([bfast],[0.5.6],[Nils Homer - Nelson Lab - UCLA, <nhomer@cs.ucla.edu>],[bfast])

AC_COPYRIGHT([See LICENSE for copyright information.])

AC_CONFIG_AUX_DIR(config)
AM_INIT_AUTOMAKE([dist-bzip2])
AC_CONFIG_SRCDIR([config.h.in])
AC_CONFIG_HEADERS([config.h])

# Sets build_cpu, build_vendor and build_os variables
AC_CANONICAL_BUILD

# Check for C and C++ compilers
AC_PROG_CC

AC_PROG_INSTALL
AC_GNU_SOURCE

# set CFLAGS and CXXFLAGS
default_CFLAGS="-Wall -g -O3 -pthread"
extended_CFLAGS=""
case "${host_cpu}-${host_os}" in 
	# Mac OS X
	i386-darwin*)
		AC_ARG_ENABLE(universal, [  --enable-universal      enable universal binary (Intel Apple only)], [is_uni=1], [is_uni=0]) 
		case $is_uni in 
			1) 
			CFLAGS="-m64"
			AC_COMPILE_IFELSE([AC_LANG_PROGRAM],
							  [extended_CFLAGS="-arch x86_64 -arch i386 -arch ppc64 -arch ppc"],
							  [extended_CFLAGS="-arch i386 -arch ppc"]);;
			0) 
			CFLAGS="-m64"
			AC_COMPILE_IFELSE([AC_LANG_PROGRAM], [extended_CFLAGS="-m64"], []);;
		esac;;
	*)
		AC_MSG_CHECKING([if gcc accepts -m64])
		CFLAGS="-m64"
		AC_COMPILE_IFELSE([AC_LANG_PROGRAM], [extended_CFLAGS="-m64"; 
						   AC_MSG_RESULT([yes])], [extended_CFLAGS="-D_FILE_OFFSET_BITS=64"; AC_MSG_RESULT([no])]);;
esac

# Define some variables
svnrev="`svn info | grep Revision:`"
if test -n "${svnrev}" ; then
	AC_REVISION("${svnrev}") 
	AC_DEFINE_UNQUOTED(SVN_REV,"${svnrev}", [Define to be the subversion revision number])
else
	AC_REVISION([$Revision: undefined $])
	AC_DEFINE(SVN_REV,"$Revision: undefined$", [Define to be the subversion revision number])
fi
					   
AC_ARG_ENABLE(intel64, [  --enable-intel64        optimize for Intel64 CPU such as Xeon and Core2],
[extended_CFLAGS="${extended_CFLAGS} -mtune=nocona"], [])
AC_ARG_ENABLE(unoptimized-sw, [  --enable-unoptimized-sw         use the unoptimized smith waterman algorithm (use for benchmarking only)],
[AC_DEFINE(UNOPTIMIZED_SMITH_WATERMAN,1,[Define 1 if we want to use the unoptimized smith waterman (use for benchmarking only)])], 
[])
CFLAGS="${default_CFLAGS} ${extended_CFLAGS}"

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
# Checks for library functions.
AC_FUNC_MALLOC 
AC_FUNC_REALLOC
AC_CHECK_LIB([m], [pow])
AC_CHECK_FUNCS([floor pow strstr strpbrk])

# Check types
AC_CHECK_TYPE([int64_t int32_t uint32_t int8_t uint8_t])
																														
# Set proprocessor SIZEOF_short_int, SIZEOF_int, SIZEOF_long_int
AC_CHECK_SIZEOF([short int])
AC_CHECK_SIZEOF([int])
AC_CHECK_SIZEOF([long int])

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([limits.h stdint.h stdlib.h string.h sys/time.h unistd.h float.h])
AC_CHECK_HEADER([argp.h], [AC_DEFINE([HAVE_ARGP_H],1,[Define to 1 if you have the <argp.h> header file])],
[AC_MSG_WARN([argp.h header file not found.
OS X is known to have this problem.
If your system does not have argp.h then bfast will
attempt to use getopt.h for command line parsing.])])


AC_CONFIG_FILES([Makefile
				 bpreprocess/Makefile
				 bmatches/Makefile
				 balign/Makefile
				 bpostprocess/Makefile
				 butil/Makefile
				 tests/Makefile
				 scripts/Makefile
				 ])

AC_OUTPUT
