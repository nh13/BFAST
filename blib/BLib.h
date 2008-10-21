#ifndef BLIB_H_
#define BLIB_H_

#include <stdio.h>
#include <stdint.h>
#include "BIndex.h"

extern int8_t DNA[5];

int32_t GetFastaHeaderLine(FILE*, BString*);
int8_t ToLower(int8_t);
int8_t ToUpper(int8_t);
int8_t GetReverseComplimentAnyCaseBase(int8_t);
int32_t ValidateBasePair(int8_t);
int32_t IsAPowerOfTwo(unsigned int);
uint32_t Log2(uint32_t);
int8_t TransformFromIUPAC(int8_t);
void CheckBIndexes(char**, int, char**, int, int, int32_t*, int32_t*, int32_t*, int32_t*, int32_t);
FILE *OpenTmpFile(BString *, BString*);
void CloseTmpFile(FILE **, BString*);
void PrintPercentCompleteShort(double);
void PrintPercentCompleteLong(double);
void PrintContigPos(FILE*, int32_t, int32_t);
int32_t UpdateRead(BString*, int);
int32_t CheckReadAgainstIndex(BIndex*, BString*);
int32_t CheckReadBase(int8_t);
uint8_t ConvertBaseToColorSpace(uint8_t, uint8_t);
uint8_t ConvertBaseAndColor(uint8_t, uint8_t);
void ConvertReadFromColorSpace(BString*, int);
void ConvertReadToColorSpace(BString*, int*);
void NormalizeRead(BString*, int*, int8_t);
void ConvertColorsToStorage(BString*, int);
int8_t ConvertColorToStorage(int8_t);
void AdjustBounds(BReferenceGenome*, int32_t*, int32_t*, int32_t*, int32_t*);
int32_t WillGenerateValidKey(BIndex*, BString*, int);
int ValidateFileNames(BString*);

#endif
