#ifndef BSORT_H_
#define BSORT_H_

#include "../blib/AlignEntries.h"

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
	AlignEntries **entriesPtr;
	int64_t low;
	int64_t mid;
	int64_t high;
	int32_t threadID;
} ThreadSortData;

void TmpFileOpen(TmpFile*, char*);
void TmpFileClose(TmpFile*);
void TmpFileInitialize(TmpFile*);
void TmpFileUpdateMetaData(TmpFile*, AlignEntries*);
void TmpFileUpdateMetaDataHelper(TmpFile*, AlignEntry*);

void MoveAllIntoTmpFile(char*, TmpFile*, char*);
void SplitEntriesAndPrint(FILE*, TmpFile*, char*, int32_t, int32_t);
void *SortAlignEntriesHelper(void*);
void *MergeAlignEntriesHelper(void*);

#endif
