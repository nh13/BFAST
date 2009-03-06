#ifndef ALIGNEDEND_H_
#define ALIGNEDEND_H_

#include <stdio.h>
#include "BLibDefinitions.h"
#include "RGBinary.h"

int AlignedEndPrint(AlignedEnd*, FILE*, int, int);
int AlignedEndRead(AlignedEnd*, FILE*, int, int);
int AlignedEndRemoveDuplicates(AlignedEnd**, int, int);
void AlignedEndQuickSort(AlignedEnd**, int, int, int, int, double*, int);
void AlignedEndMergeSort(AlignedEnd**, int, int, int, int, double*, int);
void AlignedEndCopyAtIndex(AlignedEnd*, int, AlignedEnd*, int);
int AlignedEndCompareAtIndex(AlignedEnd*, int, AlignedEnd*, int, int);
int AlignedEndGetOneRead(AlignedEnd**, FILE*);
int AlignedEndGetAll(AlignedEnd**, FILE*);
void AlignedEndCopy(AlignedEnd*, AlignedEnd*);
void AlignedEndFree(AlignedEnd*);
void AlignedEndInitialize(AlignedEnd*);
void AlignedEndCheckReference(AlignedEnd*, RGBinary*, int);
int AlignedEndGetPivot(AlignedEnd*, int, int, int);
int64_t AlignedEndGetSize(AlignedEnd*);
#endif
