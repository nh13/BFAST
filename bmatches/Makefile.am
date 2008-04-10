bin_PROGRAMS = bmatches

bmatches_SOURCES = \
	../blib/BError.c	../blib/BError.h \
	ReadInputFiles.c	ReadInputFiles.h \
	../blib/RGIndex.c	../blib/RGIndex.h \
	../blib/RGTree.c	../blib/RGTree.h \
	../blib/RGSeqPair.c	../blib/RGSeqPair.h \
	../blib/RGMatch.c	../blib/RGMatch.h \
	FindMatches.c	FindMatches.h \
	ParseInput.c	ParseInput.h \
	Definitions.h
	
##  LDADD only lists those object files which are NOT compiled from 
##  the *_SOURCES above so we need to reference libm.a here

bmatches_LDADD = -lm 
AM_CFLAGS = -Wall -ansi -pg -g -O3 

