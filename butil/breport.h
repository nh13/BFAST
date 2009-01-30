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

void PrintEntriesToBedAndWig(AlignEntries*, RGBinary*, int, int64_t, int64_t, FILE*, FILE*);
void SplitIntoTmpFilesByContig(char*, TmpFile**, int*, char*, int, int, int, int);
void SplitEntriesAndPrint(RGBinary*, FILE*, FILE*, TmpFile*, char*, int);

#endif
