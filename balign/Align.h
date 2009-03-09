#ifndef ALIGN_H_
#define ALIGN_H_
#include "../blib/BLibDefinitions.h"
#include "Definitions.h"

int AlignRGMatches(RGMatches*, RGBinary*, AlignedRead*, int32_t, int32_t, ScoringMatrix*, int32_t, int32_t, int32_t, int32_t, int32_t, int32_t);
void AlignRGMatchesOneEnd(RGMatch*, RGBinary*, AlignedEnd*, int32_t, int32_t, ScoringMatrix*, int32_t, int32_t, double*, int32_t*);
int32_t AlignExact(char*, int32_t, char*, int32_t, ScoringMatrix*, AlignedEntry*, char, int32_t, int32_t);
void AlignMismatchesOnly(char*, int32_t, char*, int32_t, ScoringMatrix*, AlignedEntry*, int32_t, char, int32_t);
void AlignFullWithBound(char*, int32_t, char*, int32_t, ScoringMatrix*, AlignedEntry*, int32_t, char, int32_t, double);
int32_t AlignRGMatchesKeepBestScore(AlignedEnd*, double);

#endif
