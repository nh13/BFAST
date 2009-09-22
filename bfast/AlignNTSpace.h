#ifndef ALIGNNTSPACE_H_
#define ALIGNNTSPACE_H_
#include "BLibDefinitions.h"
#include "Align.h"

void AlignNTSpaceUngapped(char*, int, char*, int, ScoringMatrix*, AlignedEntry*, int, uint32_t, char);
void AlignNTSpaceFull(char*, int, char*, int, ScoringMatrix*, AlignedEntry*, AlignMatrixNT**, uint32_t, char);
void AlignNTSpaceFullWithBound(char*, int, char*, int, ScoringMatrix*, AlignedEntry*, int32_t, int32_t, AlignMatrixNT**, uint32_t, char);
void FillAlignedEntryFromMatrixNTSpace(AlignedEntry*, AlignMatrixNT**, char*, int, char*, int, int, uint32_t, char, int);

#endif
