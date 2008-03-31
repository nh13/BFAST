#ifndef READINPUTFILES_H_
#define READINPUTFILES_H_

#include <stdio.h>
#include "../blib/BLibDefinitions.h"
#include "Definitions.h"

void ReadReferenceGenome(char*, RGBinary*, int, int, int, int);
char ToLower(char);
char ToUpper(char);
void InsertSequenceLetterIntoByte(unsigned char*, int, char, char);
int ReadScoringMatrix(char*, ScoringMatrix*);
#endif
