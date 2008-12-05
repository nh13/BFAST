#ifndef ALIGN_H_
#define ALIGN_H_
#include "../blib/AlignEntry.h"
#include "Definitions.h"

int Align(char*, int, char*, int, ScoringMatrix*, AlignEntry*, char, int, int);
int FillAlignEntryFromMatrix(AlignEntry*,
		AlignMatrix**,
		char*,
		int,
		char*,
		int,
		int,
		int);

#endif
