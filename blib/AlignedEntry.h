#ifndef ALIGNEDENTRY_H_
#define ALIGNEDENTRY_H_

#include <stdio.h>
#include "BLibDefinitions.h"
#include "RGBinary.h"

int AlignedEntryPrint(AlignedEntry*, FILE*, int, int);
int AlignedEntryRead(AlignedEntry*, FILE*, int, int);
int AlignedEntryRemoveDuplicates(AlignedEntry**, int, int);
void AlignedEntryQuickSort(AlignedEntry**, int, int, int, int, double*, int);
void AlignedEntryMergeSort(AlignedEntry**, int, int, int, int, double*, int);
void AlignedEntryCopyAtIndex(AlignedEntry*, int, AlignedEntry*, int);
int AlignedEntryCompareAtIndex(AlignedEntry*, int, AlignedEntry*, int, int);
int AlignedEntryGetOneRead(AlignedEntry**, FILE*);
int AlignedEntryGetAll(AlignedEntry**, FILE*);
void AlignedEntryCopy(AlignedEntry*, AlignedEntry*);
void AlignedEntryFree(AlignedEntry*);
void AlignedEntryInitialize(AlignedEntry*);
void AlignedEntryCheckReference(AlignedEntry*, RGBinary*, int);
int AlignedEntryGetPivot(AlignedEntry*, int, int, int);
int64_t AlignedEntryGetSize(AlignedEntry*);
#endif
