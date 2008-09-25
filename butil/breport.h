#ifndef BREPORT_H_
#define BREPORT_H_

#include "../blib/AlignEntries.h"

typedef struct {
	int minPos;
	int maxPos;
	int contig;
	int numEntries;
	FILE *FP;
	char *FileName;
} TmpFile;

void TmpFileOpen(TmpFile*, char*, int);
void TmpFileClose(TmpFile*);
void TmpFileInitialize(TmpFile*);

void PrintEntriesToBedAndWig(AlignEntries*, int, int64_t, int64_t, FILE*, FILE*);
int SplitIntoTmpFilesByContig(char*, TmpFile**, char*, int, int);
void SplitEntriesAndPrint(FILE*, FILE*, TmpFile*, char*, int);

#endif
