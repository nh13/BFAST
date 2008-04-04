bin_PROGRAMS = bpreprocess

bpreprocess_SOURCES = \
	ReadInputFiles.c	ReadInputFiles.h \
	GenerateTree.c	GenerateTree.h \
	GenerateIndex.c	GenerateIndex.h \
	../blib/RGIndex.c	../blib/RGIndex.h \
	../blib/RGTree.c	../blib/RGTree.h \
	../blib/BError.c	../blib/BError.h \
	ParseInput.c	ParseInput.h \
	Definitions.h
	
##  LDADD only lists those object files which are NOT compiled from 
##  the *_SOURCES above so we need to reference libm.a here

bpreprocess_LDADD = -lm 
AM_CFLAGS = -Wall -ansi -pg -g -O3


