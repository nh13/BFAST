#ifndef RGBINARY_H_
#define RGBINARY_H_

#include <stdlib.h>
#include "BLibDefinitions.h"

void RGBinaryRead(char*, RGBinary*, int, int, int, int);
void RGBinaryDelete(RGBinary*);
void RGBinaryInsertBase(unsigned char*, int, char, char);
void RGBinaryGetSequence(RGBinary*, int, int, char, int, char*, int, int*, int*);
char RGBinaryGetBase(RGBinary*, int, int);
int RGBinaryIsRepeat(RGBinary*, int, int);
int RGBinaryIsN(RGBinary*, int, int);

#endif
