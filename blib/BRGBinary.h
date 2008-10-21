#ifndef BBINARY_H_
#define BBINARY_H_

#include <stdlib.h>
#include "BLibDefinitions.h"

void BBinaryRead(char*, BBinary*, int32_t);
void BBinaryReadBinary(RGBinary*, char*);
void BBinaryWriteBinary(RGBinary*, char*);
void BBinaryDelete(RGBinary*);
void BBinaryInsertBase(uint8_t*, int32_t, int8_t);
int32_t BBinaryGetSequence(RGBinary*, int32_t, int32_t, int8_t, char**, int32_t);
void BBinaryGetReference(RGBinary*, int32_t, int32_t, int8_t, int32_t, char**, int32_t, int32_t*, int32_t*);
int8_t BBinaryGetBase(RGBinary*, int32_t, int32_t);
int32_t BBinaryIsRepeat(RGBinary*, int32_t, int32_t);
int32_t BBinaryIsBaseRepeat(int8_t);
int32_t BBinaryIsN(RGBinary*, int32_t, int32_t);
int32_t BBinaryIsBaseN(int8_t);
void BBinaryPrintInfo(char*);

#endif
