#ifndef READINPUTFILES_H_
#define READINPUTFILES_H_

#include <stdio.h>
#include <zlib.h>
#include "RGMatches.h"
#include "RGIndex.h"

int GetNextRead(FILE*, char*, char*, char*);
int GetRead(FILE*, RGMatches*, int);
int32_t GetNextReadBuffered(char**, int32_t, char*, char*, char*);
int32_t GetReadBuffered(char**,int32_t, RGMatches*, int, int32_t);
int WriteRead(FILE*, RGMatches*);
void WriteReadsToTempFile(FILE*, FILE***, char***, int, int, int, char*, int*, int32_t);
int ReadTempReadsAndOutput(gzFile, char*, gzFile, FILE*);
void ReadRGIndex(char*, RGIndex*, int);
int GetIndexFileNames(char*, int32_t, char*, char***, int32_t***);
int32_t ReadOffsets(char*, int32_t**);
int32_t GetReads(FILE*, RGMatches*, int32_t, int32_t);

#endif
