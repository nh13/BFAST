#ifndef RGBINARY_H_
#define RGBINARY_H_

#include <stdlib.h>
#include "BLibDefinitions.h"

void RGBinaryRead(char*, RGBinary*, int32_t, int32_t, int32_t, int32_t, int32_t);
void RGBinaryReadBinary(RGBinary*, char*);
void RGBinaryWriteBinary(RGBinary*, char*);
void RGBinaryDelete(RGBinary*);
void RGBinaryInsertBase(uint8_t*, int32_t, int8_t);
void RGBinaryGetSequence(RGBinary*, int32_t, int32_t, int8_t, int32_t, char*, int32_t, int32_t*, int32_t*);
int8_t RGBinaryGetBase(RGBinary*, int32_t, int32_t);
int32_t RGBinaryIsRepeat(RGBinary*, int32_t, int32_t);
int32_t RGBinaryIsBaseRepeat(int8_t);
int32_t RGBinaryIsN(RGBinary*, int32_t, int32_t);
int32_t RGBinaryIsBaseN(int8_t);

#endif
