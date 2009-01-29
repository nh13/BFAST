#ifndef ALIGN_H_
#define ALIGN_H_
#include "../blib/AlignEntry.h"
#include "Definitions.h"

int AlignRGMatches(RGMatches*, RGBinary*, AlignEntries*, int32_t, int32_t, int32_t, int32_t, ScoringMatrix*, int32_t, int32_t, int32_t, int32_t, int32_t, int32_t);
double AlignRGMatchesOneEnd(RGMatch*, RGBinary*, AlignEntry*, int32_t, int32_t, int32_t, ScoringMatrix*, int32_t, int32_t);
int32_t AlignExact(char*, int32_t, char*, int32_t, int32_t, ScoringMatrix*, AlignEntry*, char, int32_t, int32_t);
void AlignMismatchesOnly(char*, int32_t, char*, int32_t, int32_t, ScoringMatrix*, AlignEntry*, int32_t, char, int32_t);
void AlignFullWithBound(char*, int32_t, char*, int32_t, int32_t, ScoringMatrix*, AlignEntry*, int32_t, char, int32_t, double);
int32_t AlignRGMatchesKeepBestScore(AlignEntry**, int32_t*, double);

#endif
