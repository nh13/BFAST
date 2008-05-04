#ifndef RGMATCH_H_
#define RGMATCH_H_

#include <stdio.h>
#include "BLibDefinitions.h"

void RGMatchRemoveDuplicates(RGMatch*);
void RGMatchQuickSort(RGMatch*, int32_t, int32_t);
void RGMatchOutputToFile(FILE*, char*, char*, char*, RGMatch*, RGMatch*, int32_t);
int32_t RGMatchMergeFilesAndOutput(FILE**, int32_t, FILE*, int32_t);
int32_t RGMatchMergeThreadTempFilesIntoOutputTempFile(FILE**, int32_t, FILE*, int32_t);
int32_t RGMatchGetNextFromFile(FILE*, char*, char*, char*, RGMatch*, RGMatch*, int32_t);
int32_t RGMatchCompareAtIndex(RGMatch*, int32_t, RGMatch*, int32_t);
void RGMatchCopyAtIndex(RGMatch*, int32_t, RGMatch*, int32_t);
void RGMatchAllocate(RGMatch*, int32_t);
void RGMatchReallocate(RGMatch*, int32_t);

#endif

