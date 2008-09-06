#ifndef BLIB_H_
#define BLIB_H_

#include <stdio.h>
#include <stdint.h>
#include "RGIndex.h"

char ToLower(char);
char ToUpper(char);
void GetReverseComplimentAnyCase(char*, char*, int);
char GetReverseComplimentAnyCaseBase(char);
int ValidateBasePair(char);
int IsAPowerOfTwo(unsigned int);
char TransformFromIUPAC(char);
void CheckRGIndexes(char**, int, char**, int, int, int32_t*, int32_t*, int32_t*, int32_t*);
FILE *OpenTmpFile(char*, char**);
void CloseTmpFile(FILE **, char**);
void PrintPercentCompleteShort(double);
void PrintPercentCompleteLong(double);
int UpdateRead(char*, int);
int CheckReadAgainstIndex(RGIndex*, char*, int);
int CheckReadBase(char);
uint8_t ConvertToColorSpace(uint8_t, uint8_t);

#endif

