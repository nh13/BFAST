#ifndef READINPUTFILES_H_
#define READINPUTFILES_H_
#include "../blib/BLibDefinitions.h"
void ReadReferenceGenome(char*, int, RGList*, int, int, int, int, int);
void ReadGaps(char*, int**, int*, int*, int);
char ToLower(char);
int ValidateSequence(char*, int);

#endif
