#ifndef ALIGNNTSPACE_H_
#define ALIGNNTSPACE_H_
#include "../blib/AlignEntry.h"
#include "Definitions.h"

int AlignNTSpace(char*, int, char*, int, ScoringMatrix*, AlignEntry*, int);
int AlignNTSpaceFull(char*, int, char*, int, ScoringMatrix*, AlignEntry*);
int AlignNTSpaceMismatchesOnly(char*, int, char*, int, ScoringMatrix*, AlignEntry*);

#endif
