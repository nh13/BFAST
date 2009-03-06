#ifndef ALIGNNTSPACE_H_
#define ALIGNNTSPACE_H_
#include "../blib/AlignedEntry.h"
#include "Definitions.h"

void AlignNTSpaceMismatchesOnly(char*, int, char*, int, ScoringMatrix*, AlignedEntry*, char, int32_t);
void AlignNTSpaceFull(char*, int, char*, int, ScoringMatrix*, AlignedEntry*, char, int32_t);
void AlignNTSpaceFullWithBound(char*, int, char*, int, ScoringMatrix*, AlignedEntry*, char, int32_t, int32_t, int32_t);
int FillAlignedEntryFromMatrixNTSpace(AlignedEntry*, AlignMatrixNT**, char*, int, char*, int, int, int);

#endif
