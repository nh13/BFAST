#ifndef ALIGNCOLORSPACE_H_
#define ALIGNCOLORSPACE_H_
#include "../blib/AlignEntry.h"
#include "../blib/BLibDefinitions.h"
#include "Definitions.h"

int AlignColorSpace(char*, int, char*, int, ScoringMatrix*, AlignEntry*, char, int, int);
int AlignColorSpaceMismatchesOnly(char*, int, char*, int, int, ScoringMatrix*, AlignEntry*, char);
int AlignColorSpaceFull(char*, int, char*, int, int, ScoringMatrix*, AlignEntry*, char);

#endif
