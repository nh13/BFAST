#ifndef BMATCHES_H_
#define BMATCHES_H_

#include <stdio.h>
#include "BLibDefinitions.h"

int32_t BMatchesRead(FILE*, BMatches*, int32_t, int32_t);
void BMatchesPrint(FILE*, BMatches*, int32_t, int32_t);
void BMatchesRemoveDuplicates(RGMatches*, int32_t);
int32_t BMatchesMergeFilesAndOutput(FILE**, int32_t, FILE*, int32_t, int32_t, int32_t);
int32_t BMatchesMergeThreadTempFilesIntoOutputTempFile(FILE**, int32_t, FILE*, int32_t, int32_t);
int32_t BMatchesCompareAtIndex(RGMatches*, int32_t, BMatches*, int32_t);
void BMatchesAppend(RGMatches*, BMatches*);
void BMatchesCopyAtIndex(RGMatches*, int32_t, BMatches*, int32_t);
void BMatchesAllocate(RGMatches*, int32_t);
void BMatchesReallocate(RGMatches*, int32_t);
void BMatchesFree(RGMatches*);
void BMatchesInitialize(RGMatches*);
void BMatchesMirrorPairedEnd(RGMatches*, int32_t, int32_t);
void BMatchesCheck(RGMatches*);
void BMatchesFilterOutOfRange(RGMatches*, int32_t, int32_t, int32_t, int32_t, int32_t);

#endif

