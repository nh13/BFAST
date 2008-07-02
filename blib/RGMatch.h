#ifndef RGMATCH_H_
#define RGMATCH_H_

#include <stdio.h>
#include "BLibDefinitions.h"

void RGMatchRemoveDuplicates(RGMatch*, int32_t);
void RGMatchQuickSort(RGMatch*, int32_t, int32_t);
int32_t RGMatchRead(FILE*, char*, char*, char*, RGMatch*, RGMatch*, int32_t, int32_t);
void RGMatchPrint(FILE*, char*, char*, char*, RGMatch*, RGMatch*, int32_t, int32_t);
int32_t RGMatchMergeFilesAndOutput(FILE**, int32_t, FILE*, int32_t, int32_t, int32_t);
int32_t RGMatchMergeThreadTempFilesIntoOutputTempFile(FILE**, int32_t, FILE*, int32_t, int32_t);
int32_t RGMatchCompareAtIndex(RGMatch*, int32_t, RGMatch*, int32_t);
void RGMatchAppend(RGMatch*, RGMatch*, int);
void RGMatchCopyAtIndex(RGMatch*, int32_t, RGMatch*, int32_t);
void RGMatchAllocate(RGMatch*, int32_t);
void RGMatchReallocate(RGMatch*, int32_t);
void RGMatchFree(RGMatch*);
void RGMatchInitialize(RGMatch*);
void RGMatchMirrorPairedEnd(RGMatch*, RGMatch*, int);

#endif

