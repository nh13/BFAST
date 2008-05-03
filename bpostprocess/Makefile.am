bin_PROGRAMS = bpostprocess

bpostprocess_SOURCES = \
	../blib/BError.c	../blib/BError.h \
	../blib/BLib.c	../blib/BLib.h \
	../blib/RGIndex.c ../blib/RGIndex.h \
	../blib/RGBinary.c ../blib/RGBinary.h \
	../blib/AlignEntry.c	../blib/AlignEntry.h \
	PrintOutputFiles.c	PrintOutputFiles.h \
	FilterAlignments.c	FilterAlignments.h \
	ParseInput.c	ParseInput.h \
	Definitions.h
	
##  LDADD only lists those object files which are NOT compiled from 
##  the *_SOURCES above so we need to reference libm.a here

bpostprocess_LDADD = -lm 
AM_CFLAGS = -Wall -ansi -g -O3 -m64
