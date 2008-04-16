#ifndef RGREADS_H_
#define RGREADS_H_

#include <stdio.h>
#include "BLibDefinitions.h"
#include "RGMatch.h"
#include "RGIndex.h"

void RGReadsFindMatches(RGIndex*, RGBinary*, RGMatch*, char*, int*, int, int, int, int, int, int);
void RGReadsGenerateReads(char*, RGIndex*, RGReads*, char, int*, int, int, int, int, int, int);
void RGReadsGenerateMismatches(char*, int, char, int, unsigned int, unsigned int*, unsigned int*, unsigned int, int, RGReads*);
void RGReadsGenerateMismatchesHelper(char*, int, char, int, unsigned int, unsigned int*, unsigned int*, unsigned int, int, int, RGReads*, char*, int, int);
void RGReadsGenerateDeletions(char*, int, char, int, unsigned int, unsigned int*, unsigned int*, unsigned int, int, RGReads*);
void RGReadsGenerateDeletionsHelper(char*, int, char, int, unsigned int, unsigned int*, unsigned int*, unsigned int, int, int, int, int, RGReads*, char*, int, int);
void RGReadsGenerateInsertions(char*, int, char, int, unsigned int, unsigned int*, unsigned int*, unsigned int, int, RGReads*);
void RGReadsGenerateInsertionsHelper(char*, int, char, int, unsigned int, unsigned int*, unsigned int*, unsigned int, int, int, int, int, RGReads*, char*, int, int);
void RGReadsGenerateGapDeletions(char*, int, char, int, unsigned int, unsigned int*, unsigned int*, unsigned int, int, RGReads*);
void RGReadsGenerateGapDeletionsHelper(char*, int, char, int, unsigned int, unsigned int*, unsigned int*, unsigned int, int, RGReads*, char*);
void RGReadsGenerateGapInsertions(char*, int, char, int, unsigned int, unsigned int*, unsigned int*, unsigned int, int, RGReads*);
void RGReadsGenerateGapInsertionsHelper(char*, int, char, int, unsigned int, unsigned int*, unsigned int*, unsigned int, int, RGReads*, char*);

void RGReadsRemoveDuplicates(RGReads*);
void RGReadsQuickSort(RGReads*, int, int);

int RGReadsCompareAtIndex(RGReads*, int, RGReads*, int);
void RGReadsCopyAtIndex(RGReads*, int, RGReads*, int);

void RGReadsAllocate(RGReads*, int);
void RGReadsReallocate(RGReads*, int);

#endif

