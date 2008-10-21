#ifndef BRANGES_H_
#define BRANGES_H_

#include <stdio.h>
#include "BLibDefinitions.h"

void BRangesRemoveDuplicates(RGRanges*);
void BRangesCopyToRGMatch(RGRanges*, BIndex*, BMatch*);
void BRangesQuickSort(RGRanges*, int32_t, int32_t);
int32_t BRangesCompareAtIndex(RGRanges*, int32_t, BRanges*, int32_t);
void BRangesAppend(RGRanges*, BRanges*);
void BRangesCopyAtIndex(RGRanges*, int32_t, BRanges*, int32_t);
void BRangesAllocate(RGRanges*, int32_t);
void BRangesReallocate(RGRanges*, int32_t);
void BRangesFree(RGRanges*);
void BRangesInitialize(RGRanges*);

#endif

