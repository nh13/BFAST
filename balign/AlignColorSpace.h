#ifndef ALIGNCOLORSPACE_H_
#define ALIGNCOLORSPACE_H_
#include "../blib/AlignEntry.h"
#include "../blib/BLibDefinitions.h"
#include "Definitions.h"

void AlignColorSpaceMismatchesOnly(char*, int, char*, int, int, ScoringMatrix*, AlignEntry*, char, int32_t);
void AlignColorSpaceFull(char*, int, char*, int, int, ScoringMatrix*, AlignEntry*, char, int32_t);
void AlignColorSpaceFullWithBound(char*, int, char*, int, int, ScoringMatrix*, AlignEntry*, char, int32_t, int32_t, int32_t);
int FillAlignEntryFromMatrixColorSpace(AlignEntry*, AlignMatrixCS**, char*, int, char*, int, int, int, int);

#endif
