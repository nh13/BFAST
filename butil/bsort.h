#ifndef BSORT_H_
#define BSORT_H_

#include "../blib/BLibDefinitions.h"

#define BSORT_THREADED_SORT_MIN 100000

typedef struct {
	int64_t startContig;
	int64_t startPos;
	int64_t endContig;
	int64_t endPos;
	int64_t numEntries;
	FILE *FP;
	char *FileName;
} TmpFile;

typedef struct {
	AlignedRead **entriesPtr;
	int64_t low;
	int64_t mid;
	int64_t high;
	int32_t threadID;
} ThreadSortData;

void TmpFileOpen(TmpFile*, char*);
void TmpFileClose(TmpFile*);
void TmpFileInitialize(TmpFile*);
void TmpFileUpdateMetaData(TmpFile*, AlignedRead*);
void TmpFileUpdateMetaDataHelper(TmpFile*, AlignedEntry*);

void MoveAllIntoTmpFile(char*, TmpFile*, char*);
void SplitEntriesAndPrint(FILE*, TmpFile*, char*, int32_t, int32_t);
void *SortAlignedReadHelper(void*);
void *MergeAlignedReadHelper(void*);

#endif
