#ifndef READINPUTFILES_H_
#define READINPUTFILES_H_

#include <stdio.h>
#include <zlib.h>
#include "../blib/RGMatches.h"
#include "../blib/RGIndex.h"

int GetNextRead(FILE*, char*, char*, char*);
int GetRead(FILE*, RGMatches*);
int WriteRead(FILE*, RGMatches*);
void WriteReadsToTempFile(FILE*, FILE***, char***, int, int, int, char*, int*, int32_t);
int ReadTempReadsAndOutput(gzFile, char*, gzFile, FILE*);
void ReadRGIndex(char*, RGIndex*, int);
int ReadFileNames(char*, char***);
int ReadOffsets(char*, int**);

#endif
