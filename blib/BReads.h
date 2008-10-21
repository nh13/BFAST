#ifndef BREADS_H_
#define BREADS_H_

#include <stdio.h>
#include "BLibDefinitions.h"
#include "RGMatch.h"
#include "RGIndex.h"

void BReadsFindMatches(RGIndex*, BRGBinary*, BMatch*, int*, int, int, int, int, int, int, int, int, int);
void BReadsGenerateReads(char*, int, BIndex*, BReads*, char, int*, int, int, int, int, int, int, int);

void BReadsGeneratePerfectMatch(char*, int, char, int, BIndex*, BReads*);

void BReadsGenerateMismatches(char*, int, char, int, int, BIndex*, BReads*);
void BReadsGenerateMismatchesHelper(char*, int, char, int, int, char*, int, BIndex*, BReads*);

void BReadsGenerateDeletions(char*, int, char, int, int, BIndex*, BReads*);
void BReadsGenerateDeletionsHelper(char*, int, char, int, int, int, int, char*, int, BIndex*, BReads*);

void BReadsGenerateInsertions(char*, int, char, int, int, BIndex*, BReads*);
void BReadsGenerateInsertionsHelper(char*, int, char, int, int, int, int, char*, int, BIndex*, BReads*);

void BReadsGenerateGapDeletions(char*, int, char, int, int, BIndex*, BReads*);
void BReadsGenerateGapDeletionsHelper(char*, int, char, int, int, char*, BIndex*, BReads*);

void BReadsGenerateGapInsertions(char*, int, char, int, int, BIndex*, BReads*);
void BReadsGenerateGapInsertionsHelper(char*, int, char, int, int, char*, BIndex*, BReads*);

void BReadsRemoveDuplicates(RGReads*);
void BReadsQuickSort(RGReads*, int, int);

int BReadsCompareAtIndex(RGReads*, int, BReads*, int);
void BReadsCopyAtIndex(RGReads*, int, BReads*, int);

void BReadsAllocate(RGReads*, int);
void BReadsReallocate(RGReads*, int);
void BReadsFree(RGReads*);
void BReadsInitialize(RGReads*);
void BReadsAppend(RGReads*, char*, int32_t, int8_t, int32_t);
void BReadsPrint(RGReads*, BIndex*);

#endif

