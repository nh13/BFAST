#ifndef RUNALIGNER_H_
#define RUNALIGNER_H_

#include "../blib/BLibDefinitions.h"

void RunAligner(RGBinary*, char*, char*, int, int, int, char*, char*);
void RunDynamicProgramming(FILE*, RGBinary*, char*, int, int, int, FILE*);
void GetSequence(RGBinary*, int, int, char, int, int, char*);

#endif
