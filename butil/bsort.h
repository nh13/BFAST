#ifndef BSORT_H_
#define BSORT_H_

#include "../blib/AlignEntries.h"

typedef struct {
	int64_t startContig;
	int64_t startPos;
	int64_t endContig;
	int64_t endPos;
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
