#ifndef RGMATCHES_H_
#define RGMATCHES_H_

#include <stdio.h>
#include "BLibDefinitions.h"

int32_t RGMatchesRead(FILE*, RGMatches*, int32_t);
void RGMatchesPrint(FILE*, RGMatches*, int32_t);
void RGMatchesRemoveDuplicates(RGMatches*, int32_t);
int32_t RGMatchesMergeFilesAndOutput(FILE**, int32_t, FILE*, int32_t, int32_t);
int32_t RGMatchesMergeThreadTempFilesIntoOutputTempFile(FILE**, int32_t, FILE*, int32_t);
int32_t RGMatchesCompareAtIndex(RGMatches*, int32_t, RGMatches*, int32_t);
void RGMatchesAppend(RGMatches*, RGMatches*);
void RGMatchesAllocate(RGMatches*, int32_t);
void RGMatchesReallocate(RGMatches*, int32_t);
void RGMatchesFree(RGMatches*);
void RGMatchesInitialize(RGMatches*);
void RGMatchesMirrorPairedEnd(RGMatches*, RGBinary *rg, int32_t, int32_t, int32_t);
void RGMatchesCheck(RGMatches*);
void RGMatchesFilterOutOfRange(RGMatches*, int32_t, int32_t, int32_t, int32_t, int32_t);

#endif

