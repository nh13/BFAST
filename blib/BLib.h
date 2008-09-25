#ifndef BLIB_H_
#define BLIB_H_

#include <stdio.h>
#include <stdint.h>
#include "RGIndex.h"

extern char DNA[5];

int GetFastaHeaderLine(FILE*, char*);
char ToLower(char);
char ToUpper(char);
void ReverseRead(char*, char*, int);
void GetReverseComplimentAnyCase(char*, char*, int);
char GetReverseComplimentAnyCaseBase(char);
int ValidateBasePair(char);
int IsAPowerOfTwo(unsigned int);
char TransformFromIUPAC(char);
void CheckRGIndexes(char**, int, char**, int, int, int32_t*, int32_t*, int32_t*, int32_t*, int32_t);
FILE *OpenTmpFile(char*, char**);
void CloseTmpFile(FILE **, char**);
void PrintPercentCompleteShort(double);
void PrintPercentCompleteLong(double);
void PrintContigPos(FILE*, int32_t, int32_t);
int UpdateRead(char*, int);
int CheckReadAgainstIndex(RGIndex*, char*, int);
int CheckReadBase(char);
uint8_t ConvertBaseToColorSpace(uint8_t, uint8_t);
uint8_t ConvertBaseAndColor(uint8_t, uint8_t);
int ConvertReadFromColorSpace(char*, int);
void ConvertReadToColorSpace(char**, int*);
void NormalizeRead(char**, int*, char);
void ConvertColorsToStorage(char*, int);
char ConvertColorToStorage(char);
void AdjustBounds(RGBinary*, int32_t*, int32_t*, int32_t*, int32_t*);

#endif

