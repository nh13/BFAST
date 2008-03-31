#ifndef READINPUTFILES_H_
#define READINPUTFILES_H_

#include <stdio.h>
#include "../blib/BLibDefinitions.h"
#include "Definitions.h"

void ReadReferenceGenome(char*, RGBinary*, int, int, int, int);
char ToLower(char);
void InsertSequenceLetterIntoByte(char*, int, char);
int ReadScoringMatrix(char*, ScoringMatrix*);
#endif
