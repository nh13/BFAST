#  Copyright (C) 2008, University of California - Los Angeles
#
#  The version number in AC_INIT is for the bfast package and 
#  should not be used for individual binaries such as bpreprocess.
#  Developers - the version number should *NEVER* be changed except by, 
#  or with explicit permission from, Nils Homer.

AC_INIT([bfast],
        [0.1.1],
        [Nils Homer - Nelson Lab - UCLA, <nhomer@cs.ucla.edu>],
        [bfast])

AC_REVISION($Id: configure.ac,v 0.1.1 2008-07-01 00:00:00 nhomer $)
AC_COPYRIGHT(
[bfast is copyright 2008 by The University of California - Los 
Angeles.  All rights reserved.  This License is limited to, and you
may use the Software solely for, your own internal and non-commercial
use for academic and research purposes. Without limiting the foregoing,
you may not use the Software as part of, or in any way in connection 
with the production, marketing, sale or support of any commercial 
product or service or for any governmental purposes.  For commercial
or governmental use, please contact licensing@ucla.edu.  By installing 
this Software you are agreeing to the terms of the LICENSE file 
distributed with this software.  All rights for use are given to Nils Homer.])

AC_PREREQ(2.60)
AC_CONFIG_AUX_DIR(config)
AM_INIT_AUTOMAKE([dist-bzip2])
AM_CONFIG_HEADER(config.h)

# Check for C and C++ compilers
AC_PROG_CC

AC_PROG_INSTALL
AC_USE_SYSTEM_EXTENSIONS
AC_STDC_HEADERS

# Sets build_cpu, build_vendor and build_os variables
AC_CANONICAL_BUILD

AC_CHECK_LIB([m], [pow])
AC_CHECK_FUNCS([pow])
AC_C_CONST
AC_TYPE_INT32_T
AC_TYPE_INT64_T
AC_TYPE_INT8_T
AC_TYPE_UINT32_T
AC_TYPE_UINT8_T
AC_FUNC_MALLOC 
AC_FUNC_REALLOC

# Set proprocessor SIZEOF_short_int, SIZEOF_int, SIZEOF_long_int
AC_CHECK_SIZEOF([short int])
AC_CHECK_SIZEOF([int])
AC_CHECK_SIZEOF([long int])

AC_CHECK_HEADER([argp.h],
                [AC_DEFINE([HAVE_ARGP_H],1,
                           [Define to 1 if you have the <argp.h> header file])],
                [AC_MSG_WARN([argp.h header file not found.
	   OS X is known to have this problem.
           If your system does not have argp.h then bfast will
           attempt to use getopt.h for command line parsing.])])
AC_CHECK_HEADERS([limits.h],
		[AC_DEFINE([HAVE_LIMITS_H], 1,
				[Define to 1 if you have <limits.h> header file])],
		[AC_MSG_ERROR([limits.h not present])])
AC_CHECK_HEADER([stdint.h],
		[AC_DEFINE([HAVE_STDINT_H], 1,
				[Define to 1 if you have <stdint.h> header file])],
		[AC_MSG_ERROR([stdint.h not present])])
AC_CHECK_HEADER([stdlib.h],
		[AC_DEFINE([HAVE_STDLIB_H], 1,
				[Define to 1 if you have <stdlib.h> header file])],
		[AC_MSG_ERROR([stdlib.h not present])])
AC_CHECK_HEADER([string.h],
		[AC_DEFINE([HAVE_STRING_H], 1,
				[Define to 1 if you have <string.h> header file])],
		[AC_MSG_ERROR([string.h not present])])
AC_CHECK_HEADERS([sys/time.h],
		[AC_DEFINE([HAVE_SYS_TIME_H], 1,
				[Define to 1 if you have <sys/time.h> header file])],
		[AC_MSG_ERROR([sys/time.h absent. Resource usage statistic won't be reported])]) 
AC_CHECK_HEADER([unistd.h],
		[AC_DEFINE([HAVE_UNISTD_H], 1,
				[Define to 1 if you have <unistd.h> header file])],
		[AC_MSG_ERROR([unistd.h not present])])

AC_CONFIG_FILES([
   Makefile
   doc/Makefile
   bpreprocess/Makefile
   bmatches/Makefile
   balign/Makefile
   bpostprocess/Makefile
])

AC_OUTPUT
