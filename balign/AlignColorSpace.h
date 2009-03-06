#ifndef ALIGNCOLORSPACE_H_
#define ALIGNCOLORSPACE_H_
#include "../blib/AlignedEntry.h"
#include "../blib/BLibDefinitions.h"
#include "Definitions.h"

void AlignColorSpaceMismatchesOnly(char*, int, char*, int, int, ScoringMatrix*, AlignedEntry*, char, int32_t);
void AlignColorSpaceFull(char*, int, char*, int, int, ScoringMatrix*, AlignedEntry*, char, int32_t);
void AlignColorSpaceFullWithBound(char*, int, char*, int, int, ScoringMatrix*, AlignedEntry*, char, int32_t, int32_t, int32_t);
int FillAlignedEntryFromMatrixColorSpace(AlignedEntry*, AlignMatrixCS**, char*, int, char*, int, int, int, int);

#endif
