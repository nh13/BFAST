#ifndef BSORT_H_
#define BSORT_H_

#include <zlib.h>
#include "../bfast/BLibDefinitions.h"

#define BSORT_THREADED_SORT_MIN 100000

typedef struct {
	int64_t startContig;
	int64_t startPos;
	int64_t endContig;
	int64_t endPos;
	int64_t numEntries;
	gzFile FP;
	char *FileName;
} TmpGZFile;

typedef struct {
	AlignedRead **entriesPtr;
	int64_t low;
	int64_t mid;
	int64_t high;
	int32_t threadID;
} ThreadSortData;

void TmpGZFileOpen(TmpGZFile*, char*);
void TmpGZFileClose(TmpGZFile*);
void TmpGZFileInitialize(TmpGZFile*);
void TmpGZFileUpdateMetaData(TmpGZFile*, AlignedRead*);
void TmpGZFileUpdateMetaDataHelper(TmpGZFile*, AlignedEntry*);

void MoveAllIntoTmpGZFile(char*, TmpGZFile*, char*);
void SplitEntriesAndPrint(gzFile, TmpGZFile*, char*, int32_t, int32_t);
void *SortAlignedReadHelper(void*);
void *MergeAlignedReadHelper(void*);

#endif
