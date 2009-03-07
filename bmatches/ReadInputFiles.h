#ifndef READINPUTFILES_H_
#define READINPUTFILES_H_

#include <stdio.h>
#include "../blib/RGMatches.h"
#include "../blib/RGIndex.h"

int GetNextRead(FILE*, char*, char*, char*);
int GetRead(FILE*, RGMatches*);
int WriteRead(FILE*, RGMatches*);
void WriteReadsToTempFile(FILE*, FILE*, FILE***, char***, int, int, int, char*, int*, int*);
int ReadTempReadsAndOutput(FILE*, FILE*, FILE*, int);
void ReadRGIndex(char*, RGIndex*, int);
int ReadFileNames(char*, char***);
int ReadOffsets(char*, int**);

#endif
