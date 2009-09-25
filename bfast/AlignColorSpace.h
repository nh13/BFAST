#ifndef ALIGNCOLORSPACE_H_
#define ALIGNCOLORSPACE_H_
#include "AlignedEntry.h"
#include "BLibDefinitions.h"

int32_t AlignColorSpaceUngapped(char*, char*, int, char*, int, int, ScoringMatrix*, AlignedEntry*, int, int32_t, char);
void AlignColorSpaceFull(char*, int, char*, int, ScoringMatrix*, AlignedEntry*, AlignMatrix*, int32_t, char);
void AlignColorSpaceGappedBounded(char*, int, char*, int, ScoringMatrix*, AlignedEntry*, AlignMatrix*, int32_t, char, int32_t, int32_t);
void FillAlignedEntryFromMatrixColorSpace(AlignedEntry*, AlignMatrix*, char*, int, char*, int, int, int32_t, char, int, int);

#endif
