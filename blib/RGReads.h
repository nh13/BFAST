#ifndef RGREADS_H_
#define RGREADS_H_

#include <stdio.h>
#include "BLibDefinitions.h"
#include "RGMatch.h"
#include "RGIndex.h"

void RGReadsFindMatches(RGIndex*, RGBinary*, RGMatch*, int*, int, int, int, int, int, int, int, int);
void RGReadsGenerateReads(char*, int, RGIndex*, RGReads*, char, int*, int, int, int, int, int, int, int);
void RGReadsGeneratePerfectMatch(char*, int, char, int, int32_t, int32_t*, int32_t*, int64_t, RGReads*);
void RGReadsGenerateMismatches(char*, int, char, int, int32_t, int32_t*, int32_t*, int64_t, int, RGReads*);
void RGReadsGenerateMismatchesHelper(char*, int, char, int, int32_t, int32_t*, int32_t*, int64_t, int, int, RGReads*, char*, int, int);
void RGReadsGenerateDeletions(char*, int, char, int, int32_t, int32_t*, int32_t*, int64_t, int, RGReads*);
void RGReadsGenerateDeletionsHelper(char*, int, char, int, int32_t, int32_t*, int32_t*, int64_t, int, int, int, int, RGReads*, char*, int, int);
void RGReadsGenerateInsertions(char*, int, char, int, int32_t, int32_t*, int32_t*, int64_t, int, RGReads*);
void RGReadsGenerateInsertionsHelper(char*, int, char, int, int32_t, int32_t*, int32_t*, int64_t, int, int, int, int, RGReads*, char*, int, int);
void RGReadsGenerateGapDeletions(char*, int, char, int, int32_t, int32_t*, int32_t*, int64_t, int, RGReads*);
void RGReadsGenerateGapDeletionsHelper(char*, int, char, int, int32_t, int32_t*, int32_t*, int64_t, int, RGReads*, char*);
void RGReadsGenerateGapInsertions(char*, int, char, int, int32_t, int32_t*, int32_t*, int64_t, int, RGReads*);
void RGReadsGenerateGapInsertionsHelper(char*, int, char, int, int32_t, int32_t*, int32_t*, int64_t, int, RGReads*, char*);

void RGReadsRemoveDuplicates(RGReads*);
void RGReadsQuickSort(RGReads*, int, int);

int RGReadsCompareAtIndex(RGReads*, int, RGReads*, int);
void RGReadsCopyAtIndex(RGReads*, int, RGReads*, int);

void RGReadsAllocate(RGReads*, int);
void RGReadsReallocate(RGReads*, int);
void RGReadsFree(RGReads*);
void RGReadsInitialize(RGReads*);
void RGReadsAppend(RGReads*, char*, int32_t, int8_t, int32_t);
void RGReadsPrint(RGReads*, RGIndex*);

#endif

