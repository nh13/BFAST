#ifndef ALIGNCOLORSPACE_H_
#define ALIGNCOLORSPACE_H_
#include "AlignedEntry.h"
#include "BLibDefinitions.h"

void AlignColorSpaceUngapped(char*, int, char*, int, ScoringMatrix*, AlignedEntry*, int, uint32_t, char);
void AlignColorSpaceFull(char*, int, char*, int, ScoringMatrix*, AlignedEntry*, AlignMatrixCS**, uint32_t, char);
void AlignColorSpaceFullWithBound(char*, int, char*, int, ScoringMatrix*, AlignedEntry*, int32_t, int32_t, AlignMatrixCS**, uint32_t, char);
void FillAlignedEntryFromMatrixColorSpace(AlignedEntry*, AlignMatrixCS**, char*, int, char*, int, int, uint32_t, char, int);

#endif
