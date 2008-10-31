#ifndef RGRANGES_H_
#define RGRANGES_H_

#include <stdio.h>
#include "BLibDefinitions.h"

int RGRangesRemoveDuplicates(RGRanges*);
void RGRangesCopyToRGMatch(RGRanges*, RGIndex*, RGMatch*);
void RGRangesQuickSort(RGRanges*, int32_t, int32_t);
int32_t RGRangesCompareAtIndex(RGRanges*, int32_t, RGRanges*, int32_t);
void RGRangesAppend(RGRanges*, RGRanges*);
void RGRangesCopyAtIndex(RGRanges*, int32_t, RGRanges*, int32_t);
void RGRangesAllocate(RGRanges*, int32_t);
void RGRangesReallocate(RGRanges*, int32_t);
void RGRangesFree(RGRanges*);
void RGRangesInitialize(RGRanges*);

#endif

