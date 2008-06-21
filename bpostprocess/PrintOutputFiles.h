#ifndef PRINTOUTPUTFILES_H_
#define PRINTOUTPUTFILES_H_

#include <stdio.h>
#include "../blib/AlignEntry.h"
#include "Definitions.h"

int PrintAlignEntriesToTempFilesByChr(FILE*,
		RGFiles*,
		int,
		int,
		int,
		int,
		int,
		int,
		int,
		char*);

int PrintAlignEntriesToTempFilesWithinChr(FILE*,
		int,
		int,
		int,
		int,
		int,
		int,
		int,
		char*,
		ChrFiles*);

void PrintAlignEntries(ChrFiles*,
		int,
		int,
		char*,
		char*,
		int);

void PrintSortedAlignEntriesToWig(AlignEntry*,
		int,
		int,
		FILE*);

void PrintSortedAlignEntriesToBed(AlignEntry*,
		int,
		int,
		FILE**,
		int,
		int);
#endif
