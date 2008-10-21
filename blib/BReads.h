#ifndef BREADS_H_
#define BREADS_H_

#include <stdio.h>
#include "BLibDefinitions.h"
#include "RGMatch.h"
#include "RGIndex.h"

void BReadsFindMatches(RGIndex*, BRGBinary*, BMatch*, int*, int, int, int, int, int, int, int, int, int);
void BReadsGenerateReads(Bstring*, BIndex*, BReads*, char, int*, int, int, int, int, int, int, int);

void BReadsGeneratePerfectMatch(Bstring*, char, int, BIndex*, BReads*);

void BReadsGenerateMismatches(Bstring*, char, int, int, BIndex*, BReads*);
void BReadsGenerateMismatchesHelper(Bstring*, char, int, int, BString*, int, BIndex*, BReads*);

void BReadsGenerateDeletions(Bstring*, char, int, int, BIndex*, BReads*);
void BReadsGenerateDeletionsHelper(Bstring*, char, int, int, int, int, BString*, int, BIndex*, BReads*);

void BReadsGenerateInsertions(Bstring*, char, int, int, BIndex*, BReads*);
void BReadsGenerateInsertionsHelper(Bstring*, char, int, int, int, int, BString*, int, BIndex*, BReads*);

void BReadsGenerateGapDeletions(Bstring*, char, int, int, BIndex*, BReads*);
void BReadsGenerateGapDeletionsHelper(Bstring*, char, int, int, BString*, BIndex*, BReads*);

void BReadsGenerateGapInsertions(Bstring*, char, int, int, BIndex*, BReads*);
void BReadsGenerateGapInsertionsHelper(Bstring*, char, int, int, BString*, BIndex*, BReads*);

void BReadsRemoveDuplicates(RGReads*);
void BReadsQuickSort(RGReads*, int, int);

int BReadsCompareAtIndex(RGReads*, int, BReads*, int);
void BReadsCopyAtIndex(RGReads*, int, BReads*, int);

void BReadsAllocate(RGReads*, int);
void BReadsReallocate(RGReads*, int);
void BReadsFree(RGReads*);
void BReadsInitialize(RGReads*);
void BReadsAppend(RGReads*, BString*, int32_t, int8_t, int32_t);
void BReadsPrint(RGReads*, BIndex*);

#endif

