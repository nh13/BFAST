#ifndef ALIGNNTSPACE_H_
#define ALIGNNTSPACE_H_
#include "../blib/AlignEntry.h"
#include "Definitions.h"

void AlignNTSpaceMismatchesOnly(char*, int, char*, int, ScoringMatrix*, AlignEntry*, char, int32_t);
void AlignNTSpaceFull(char*, int, char*, int, ScoringMatrix*, AlignEntry*, char, int32_t);
void AlignNTSpaceFullWithBound(char*, int, char*, int, ScoringMatrix*, AlignEntry*, char, int32_t, int32_t, int32_t);
int FillAlignEntryFromMatrixNTSpace(AlignEntry*, AlignMatrixNT**, char*, int, char*, int, int, int, int);

#endif
