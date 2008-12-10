#ifndef BSORT_H_
#define BSORT_H_

#include "../blib/AlignEntries.h"

typedef struct {
	int32_t startContig;
	int32_t startPos;
	int32_t endContig;
	int32_t endPos;
	int32_t memory; /* in MB */
	int64_t numEntries;
	FILE *FP;
	char *FileName;
} TmpFile;

void TmpFileOpen(TmpFile*, char*);
void TmpFileClose(TmpFile*);
void TmpFileInitialize(TmpFile*);
void TmpFileUpdateMetaData(TmpFile*, AlignEntries*);
void TmpFileUpdateMetaDataHelper(TmpFile*, AlignEntry*);

void MoveAllIntoTmpFile(char*, TmpFile*, char*);
void SplitEntriesAndPrint(FILE*, TmpFile*, char*, int32_t);

#endif
