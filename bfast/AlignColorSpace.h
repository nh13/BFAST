#ifndef ALIGNCOLORSPACE_H_
#define ALIGNCOLORSPACE_H_
#include "AlignedEntry.h"
#include "BLibDefinitions.h"

void AlignColorSpaceUngapped(char*, int, char*, int, ScoringMatrix*, AlignedEntry*, char, int32_t);
void AlignColorSpaceFull(char*, int, char*, int, ScoringMatrix*, AlignedEntry*, char, int32_t);
void AlignColorSpaceFullWithBound(char*, int, char*, int, ScoringMatrix*, AlignedEntry*, char, int32_t, int32_t, int32_t);
int FillAlignedEntryFromMatrixColorSpace(AlignedEntry*, AlignMatrixCS**, char*, int, char*, int, int, int);

#endif
