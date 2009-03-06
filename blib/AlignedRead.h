#ifndef ALIGNEDREAD_H_
#define ALIGNEDREAD_H_

#include <stdio.h>
#include "BLibDefinitions.h"
#include "RGBinary.h"

void AlignedReadPrint(AlignedRead*, FILE*, int32_t);
int AlignedReadRead(AlignedRead*, FILE*, int32_t);
void AlignedReadRemoveDuplicates(AlignedRead*, int32_t);
void AlignedReadQuickSort(AlignedRead*, int32_t, int32_t);
void AlignedReadMergeSort(AlignedRead*, int32_t, int32_t);
void AlignedReadReallocate(AlignedRead*, int32_t);
void AlignedReadAllocate(AlignedRead*, char*, int32_t, int32_t);
void AlignedReadFree(AlignedRead*);
void AlignedReadInitialize(AlignedRead*);
void AlignedReadCopy(AlignedRead*, AlignedRead*);
void AlignedReadMergeSortAll(AlignedRead**, int64_t, int64_t);
void AlignedReadMergeAll(AlignedRead**, int64_t, int64_t, int64_t);
int32_t AlignedReadCompareAll(AlignedRead*, AlignedRead*);
#endif

