bin_PROGRAMS = bmatches

bmatches_SOURCES = \
	ReadInputFiles.c	ReadInputFiles.h \
	../blib/BError.c	../blib/BError.h \
	../blib/BLib.c	../blib/BLib.h \
	../blib/RGIndex.c	../blib/RGIndex.h \
	../blib/RGBinary.c	../blib/RGBinary.h \
	../blib/RGReads.c	../blib/RGReads.h \
	../blib/RGMatch.c	../blib/RGMatch.h \
	FindMatches.c	FindMatches.h \
	ParseInput.c	ParseInput.h \
	Definitions.h
	
##  LDADD only lists those object files which are NOT compiled from 
##  the *_SOURCES above so we need to reference libm.a here

bmatches_LDADD = -lm 
AM_CFLAGS = -Wall -ansi -g -O3 -pthread -m64
