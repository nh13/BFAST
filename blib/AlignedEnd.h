#ifndef ALIGNEDEND_H_
#define ALIGNEDEND_H_

#include <stdio.h>

#include "BLibDefinitions.h"
#include "RGBinary.h"

int32_t AlignedEndPrint(AlignedEnd*, FILE*, int32_t, int32_t);
int32_t AlignedEndRead(AlignedEnd*, FILE*, int32_t, int32_t);
int32_t AlignedEndRemoveDuplicates(AlignedEnd*, int32_t);
void AlignedEndQuickSort(AlignedEnd*, int32_t, int32_t);
void AlignedEndMergeSort(AlignedEnd*, int32_t, int32_t);
int32_t AlignedEndCompare(AlignedEnd*, AlignedEnd*, int32_t);
void AlignedEndCopyAtIndex(AlignedEnd*, int32_t, AlignedEnd*, int32_t);
void AlignedEndCopy(AlignedEnd*, AlignedEnd*);
void AlignedEndAllocate(AlignedEnd*, char*, char*, int32_t);
void AlignedEndReallocate(AlignedEnd*, int32_t);
void AlignedEndFree(AlignedEnd*);
void AlignedEndInitialize(AlignedEnd*);

#endif
