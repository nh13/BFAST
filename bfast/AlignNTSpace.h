#ifndef ALIGNNTSPACE_H_
#define ALIGNNTSPACE_H_
#include "BLibDefinitions.h"
#include "AlignMatrix.h"
#include "Align.h"

void AlignNTSpaceUngapped(char*, char*, int, char*, int, int, ScoringMatrix*, AlignedEntry*, int, int32_t, char);
void AlignNTSpaceGappedBounded(char*, int, char*, int, ScoringMatrix*, AlignedEntry*, AlignMatrix*, int32_t, char, int32_t, int32_t);
void AlignNTSpaceGappedConstrained(char*, char*, int, char*, int, ScoringMatrix*, AlignedEntry*, AlignMatrix*, int, int32_t, int32_t, char);
void AlignNTSpaceGappedRun(char*, int, char*, int, ScoringMatrix*, AlignMatrix*, int32_t, int32_t, int32_t, int32_t);
void AlignNTSpaceRecoverAlignmentFromMatrix(AlignedEntry*, AlignMatrix*, char*, int, char*, int, int, int32_t, char, int);
#endif
