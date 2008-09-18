#ifndef ALIGN_H_
#define ALIGN_H_
#include "../blib/AlignEntry.h"
#include "Definitions.h"

#define ALIGN_DEBUG_ON 0

int Align(char*, int, char*, int, ScoringMatrix*, AlignEntry*, char, int);
double GetNTScore(char, char, ScoringMatrix*);
double GetColorScore(uint8_t, uint8_t, ScoringMatrix*);
double GetScoreFromMatrix(char, char, ScoringMatrix*);
int FillAlignEntryFromMatrix(AlignEntry*,
		AlignMatrix**,
		char*,
		int,
		char*,
		int,
		int);

#endif
