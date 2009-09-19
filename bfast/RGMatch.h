#ifndef RGMATCH_H_
#define RGMATCH_H_

#include <stdio.h>
#include <zlib.h>
#include "BLibDefinitions.h"

int32_t RGMatchRead(gzFile, RGMatch*);
int32_t RGMatchReadText(FILE*, RGMatch*);
void RGMatchPrint(gzFile, RGMatch*);
void RGMatchPrintText(FILE*, RGMatch*);
void RGMatchPrintFastq(FILE*, char*, RGMatch*);
void RGMatchRemoveDuplicates(RGMatch*, int32_t);
void RGMatchQuickSort(RGMatch*, int32_t, int32_t);
int32_t RGMatchCompareAtIndex(RGMatch*, int32_t, RGMatch*, int32_t);
void RGMatchAppend(RGMatch*, RGMatch*);
void RGMatchCopyAtIndex(RGMatch*, int32_t, RGMatch*, int32_t);
void RGMatchAllocate(RGMatch*, int32_t);
void RGMatchReallocate(RGMatch*, int32_t);
void RGMatchClearMatches(RGMatch*);
void RGMatchFree(RGMatch*);
void RGMatchInitialize(RGMatch*);
void RGMatchCheck(RGMatch*);
void RGMatchFilterOutOfRange(RGMatch*, int32_t);

#endif

