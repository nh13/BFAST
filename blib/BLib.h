#ifndef BLIB_H_
#define BLIB_H_

#include <stdio.h>
#include <stdint.h>
#include "RGIndex.h"
#include "BLibDefinitions.h"

extern char DNA[5];

int GetFastaHeaderLine(FILE*, char*);
char ToLower(char);
void ToLowerRead(char*, int);
char ToUpper(char);
void ReverseRead(char*, char*, int);
void GetReverseComplimentAnyCase(char*, char*, int);
char GetReverseComplimentAnyCaseBase(char);
int ValidateBasePair(char);
int IsAPowerOfTwo(unsigned int);
uint32_t Log2(uint32_t);
char TransformFromIUPAC(char);
void CheckRGIndexes(char**, int, char**, int, int32_t*, int32_t*, int32_t*, int32_t*, int32_t);
FILE *OpenTmpFile(char*, char**);
void CloseTmpFile(FILE **, char**);
void PrintPercentCompleteShort(double);
void PrintPercentCompleteLong(double);
int PrintContigPos(FILE*, int32_t, int32_t);
int UpdateRead(char*, int);
int CheckReadAgainstIndex(RGIndex*, char*, int);
int CheckReadBase(char);
int ConvertBaseToColorSpace(char, char, char*);
int ConvertBaseAndColor(char, char, char*);
int ConvertReadFromColorSpace(char*, int);
void ConvertReadToColorSpace(char**, int*);
void NormalizeRead(char**, int*, char);
void NormalizeColorSpaceRead(char *read, int, char);
void ConvertColorsToStorage(char*, int);
char ConvertColorToStorage(char);
void ConvertColorsFromStorage(char*, int);
char ConvertColorFromStorage(char);
char ConvertIntColorToCharColor(char);
void AdjustBounds(RGBinary*, int32_t*, int32_t*, int32_t*, int32_t*);
int WillGenerateValidKey(RGIndex*, char*, int);
int ValidateFileName(char*);
void StringCopyAndReallocate(char**, const char*);
int StringTrimWhiteSpace(char*);
int IsWhiteSpace(char);
void CheckPackageCompatibility(char*, int);
void KnuthMorrisPrattCreateTable(char*, int, int*);
int32_t KnuthMorrisPratt(char*, int, char*, int);
int32_t NaiveSubsequence(char*, int, char*, int);
int CompareContigPos(int32_t, int32_t, int32_t, int32_t);
int WithinRangeContigPos(int32_t, int32_t, int32_t, int32_t, int32_t, int32_t);
char *StrStrGetLast(char*, const char*);
void ParseRange(Range*, char*);
int32_t CheckRange(Range*, int32_t, int32_t);
int32_t CheckRangeWithinRange(Range*, Range*);
void RangeCopy(Range*, Range*);
int GetNumMismatchesInAlignedEntry(AlignedEntry *a);
int GetNumColorErrorsInAlignedEntry(AlignedEntry *a, int space);
double AddLog10(double, double);

#endif
