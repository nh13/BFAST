#ifndef READINPUTFILES_H_
#define READINPUTFILES_H_

#include <stdio.h>
#include "../blib/RGMatches.h"
#include "../blib/RGIndex.h"

int GetNextRead(FILE*, RGMatches*, int);
int WriteRead(FILE*, RGMatches*, int);
void WriteReadsToTempFile(FILE*, FILE***, char***, int, int, int, int, char*, int, int*, int*);
int ReadTempReadsAndOutput(FILE*, FILE*, FILE*, int, int);
void ReadRGIndex(char*, RGIndex*, int);
int ReadFileNames(char*, char***);
int ReadOffsets(char*, int**);

#endif
