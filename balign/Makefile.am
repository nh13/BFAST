bin_PROGRAMS = balign

balign_SOURCES = \
	RunAligner.c	RunAligner.h \
	ParseInput.c	ParseInput.h \
	ReadInputFiles.c	ReadInputFiles.h \
	../blib/BError.c	../blib/BError.h \
	../blib/BLib.c	../blib/BLib.h \
	../blib/RGMatch.c	../blib/RGMatch.h \
	../blib/RGBinary.c	../blib/RGBinary.h \
	../blib/RGIndex.c	../blib/RGIndex.h \
	../blib/AlignEntry.c	../blib/AlignEntry.h \
	Align.c		Align.h \
	Definitions.h
	
##  LDADD only lists those object files which are NOT compiled from 
##  the *_SOURCES above so we need to reference libm.a here

balign_LDADD = -lm 
AM_CFLAGS = -Wall -ansi -g -g -O3 -pthread -arch x86_64

