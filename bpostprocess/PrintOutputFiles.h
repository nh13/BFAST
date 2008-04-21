#ifndef PRINTOUTPUTFILES_H_
#define PRINTOUTPUTFILES_H_

#include <stdio.h>
#include "../blib/AlignEntry.h"
#include "Definitions.h"

void PrintAlignEntriesToTempFiles(FILE*,
		int,
		int,
		int,
		int,
		int,
		int,
		int,
		ChrFiles*);

void PrintAlignEntries(ChrFiles*,
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
		int);
#endif
