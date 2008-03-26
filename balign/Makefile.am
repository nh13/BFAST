bin_PROGRAMS = balign

balign_SOURCES = \
	ParseInput.c	ParseInput.h \
	BError.c	BError.h \
	Definitions.h
	
##  LDADD only lists those object files which are NOT compiled from 
##  the *_SOURCES above so we need to reference libm.a here

balign_LDADD = -lm 
AM_CFLAGS = -Wall -ansi -pg -g -O3

