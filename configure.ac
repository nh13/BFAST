#  Copyright (C) 2008, University of California - Los Angeles
#
#  The version number in AC_INIT is for the Blatter package and 
#  should not be used for individual binaries such as bpreprocess.
#  Developers - the version number should *NEVER* be changed except by, 
#  or with explicit permission from, Nils Homer.

AC_INIT([Blatter],
        [0.1.1],
        [Nils Homer - Nelson Lab - UCLA, nhomer@cs.ucla.edu],
        [Blatter])

AC_REVISION($Id: configure.ac,v 1.1.1.1 2007-09-11 08:36:59 nhomer Exp $)
AC_COPYRIGHT(
[Blatter is copyright 2008 by The University of California - Los 
Angeles.  All rights reserved.  This License is limited to, and you
may use the Software solely for, your own internal and non-commercial
use for academic and research purposes. Without limiting the foregoing,
you may not use the Software as part of, or in any way in connection 
with the production, marketing, sale or support of any commercial 
product or service or for any governmental purposes.  For commercial
or governmental use, please contact licensing@tgen.org.  By installing 
this Software you are agreeing to the terms of the LICENSE file 
distributed with this software.  All rights for use are given to Nils Homer.])

AC_CONFIG_AUX_DIR(config)
AM_INIT_AUTOMAKE([dist-bzip2])
AM_CONFIG_HEADER(config.h)

# Check for C and C++ compilers
AC_PROG_CC
AC_PROG_CXX

AC_PROG_INSTALL
AC_AIX
AC_ISC_POSIX
AC_MINIX
AC_STDC_HEADERS
AC_PROG_RANLIB

# Sets build_cpu, build_vendor and build_os variables
AC_CANONICAL_BUILD

# 4 extra macros added based on suggestions from running "autoscan" 
AC_CHECK_FUNCS([pow])
AC_CHECK_FUNCS([sqrt])
AC_C_CONST
AC_FUNC_MALLOC 
AC_FUNC_STRTOD
AC_HEADER_STDBOOL
AC_CHECK_FUNCS([strchr])
AC_CHECK_FUNCS([strstr])

# Set proprocessor SIZEOF_short_int, SIZEOF_int, SIZEOF_long_int
AC_CHECK_SIZEOF([short int])
AC_CHECK_SIZEOF([int])
AC_CHECK_SIZEOF([long int])

AC_CHECK_HEADER([argp.h],
                [AC_DEFINE([HAVE_ARGP_H],1,
                           [Define to 1 if you have the <argp.h> header file])],
                [AC_MSG_WARN([argp.h header file not found.
	   OS X is known to have this problem.
           If your system does not have argp.h then Blatter will
           attempt to use getopt.h for command line parsing.])])

AC_CHECK_HEADER([netinet/in.h],
		[AC_DEFINE([HAVE_NETINET_H], 1,
				[Define to 1 if htons/ntohs expected in <netinet/in.h>])],
		[AC_MSG_WARN([netinet/in.h not present])])
AC_CHECK_HEADER([limits.h],
		[AC_DEFINE([HAVE_LIMITS_H], 1,
				[Define to 1 if you have <limits.h> header file])],
		[AC_MSG_WARN([limits.h not present])])
AC_CHECK_HEADER([libintl.h],
		[AC_DEFINE([HAVE_LIBINTL_H], 1,
				[Define to 1 if you have libintl.h header file])],
		[AC_MSG_WARN([libintl.h not present])])
AC_CHECK_HEADER([sys/resource.h],
		[AC_DEFINE([HAVE_RESOURCE_H], 1, 
				[Define to 1 if you have <sys/resource.h> header file])],
		[AC_MSG_WARN([resource.h absent. Resource usage statistics won't be reported])])
AC_CHECK_HEADERS([sys/time.h],
		[AC_DEFINE([HAVE_SYS_TIME_H], 1,
				[Define to 1 if you have <sys/time.h> header file])],
		[AC_MSG_WARN([sys/time.h absent. Resource usage statistic won't be reported])]) 
AC_CHECK_HEADERS([sys/types.h],
		[AC_DEFINE([HAVE_SYS_TYPES_H], 1,
				[Define to 1 if you have <sys/type.h> header file])],
		[AC_MSG_WARN([sys/type.h absent.])]) 
AC_CHECK_HEADERS([sys/socket.h],
		[AC_DEFINE([HAVE_SYS_SOCKET_H], 1,
				[Define to 1 if you have <sys/socket.h> header file])],
		[AC_MSG_WARN([sys/socket.h absent.])])

AC_CONFIG_FILES([
   Makefile
   doc/Makefile
   bpreprocess/Makefile
   bmatches/Makefile
   balign/Makefile
])

AC_OUTPUT
