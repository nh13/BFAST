bin_PROGRAMS = bpreprocess

bpreprocess_SOURCES = \
	GenerateIndex.c	GenerateIndex.h \
	../blib/BLib.c	../blib/BLib.h \
	../blib/RGBinary.c	../blib/RGBinary.h \
	../blib/RGIndexLayout.c	../blib/RGIndexLayout.h \
	../blib/RGIndex.c	../blib/RGIndex.h \
	../blib/RGMatch.c	../blib/RGMatch.h \
	../blib/BError.c	../blib/BError.h \
	ParseInput.c	ParseInput.h \
	Definitions.h
	
##  LDADD only lists those object files which are NOT compiled from 
##  the *_SOURCES above so we need to reference libm.a here

AM_CFLAGS = -Wall -ansi -g -O0 -pthread -m64
