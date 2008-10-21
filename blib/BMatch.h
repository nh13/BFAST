#ifndef BMATCH_H_
#define BMATCH_H_

#include <stdio.h>
#include "BLibDefinitions.h"

int32_t BMatchRead(FILE*, BMatch*, int32_t);
void BMatchPrint(FILE*, BMatch*, int32_t);
void BMatchRemoveDuplicates(RGMatch*, int32_t);
void BMatchQuickSort(RGMatch*, int32_t, int32_t);
int32_t BMatchCompareAtIndex(RGMatch*, int32_t, BMatch*, int32_t);
void BMatchAppend(RGMatch*, BMatch*);
void BMatchCopyAtIndex(RGMatch*, int32_t, BMatch*, int32_t);
void BMatchAllocate(RGMatch*, int32_t);
void BMatchReallocate(RGMatch*, int32_t);
void BMatchClearMatches(RGMatch*);
void BMatchFree(RGMatch*);
void BMatchInitialize(RGMatch*);
void BMatchCheck(RGMatch*);
void BMatchFilterOutOfRange(RGMatch*, int32_t, int32_t, int32_t, int32_t, int32_t);

#endif

