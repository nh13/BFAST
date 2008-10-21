#ifndef BINDEX_H_
#define BINDEX_H_

#include <stdio.h>
#include "BRGBinary.h"
#include "BRanges.h"
#include "BLibDefinitions.h"

void BIndexCreate(BIndex*, BRGBinary*, BIndexLayout*, int32_t, int32_t, int32_t, int32_t, int32_t, int32_t, BIndexExons*, int32_t, int32_t, int32_t, int32_t, char*);
void BIndexCreateHelper(BIndex*, BRGBinary*, int32_t, int32_t, int32_t, int32_t, int32_t, int32_t);
void BIndexCreateHash(BIndex*, BRGBinary*);
void BIndexSort(BIndex*, BRGBinary*, int32_t, char*);
void *BIndexMergeSort(void*);
void BIndexMergeSortHelper(BIndex*, BRGBinary*, int64_t, int64_t, int32_t, double*, int64_t, int64_t, int64_t, char*);
void *BIndexMerge(void*);
void BIndexMergeHelper(BIndex*, BRGBinary*, int64_t, int64_t, int64_t, int64_t, char*);
void BIndexMergeHelperInMemoryContig_8(BIndex*, BRGBinary*, int64_t, int64_t, int64_t);
void BIndexMergeHelperInMemoryContig_32(BIndex*, BRGBinary*, int64_t, int64_t, int64_t);
void BIndexMergeHelperFromDiskContig_8(BIndex*, BRGBinary*, int64_t, int64_t, int64_t, char*);
void BIndexMergeHelperFromDiskContig_32(BIndex*, BRGBinary*, int64_t, int64_t, int64_t, char*);

void BIndexDelete(BIndex*);
double BIndexGetSize(BIndex*, int32_t);
void BIndexPrint(FILE*, BIndex*, int32_t);
void BIndexRead(FILE*, BIndex*, int32_t);
void BIndexPrintInfo(char*);
void BIndexPrintHeader(FILE*, BIndex*, int32_t);
void BIndexReadHeader(FILE*, BIndex*, int32_t);
void BIndexGetRanges(BIndex*, BRGBinary*, char*, int32_t, int8_t, int32_t, int32_t, BRanges*);
int64_t BIndexGetIndex(BIndex*, BRGBinary*, int64_t, int64_t, char*, int64_t*, int64_t*);
void BIndexSwapAt(BIndex*, int64_t, int64_t);
int64_t BIndexGetPivot(BIndex*, BRGBinary*, int64_t, int64_t);
int32_t BIndexCompareContigPos(BIndex*, BRGBinary*, uint32_t, uint32_t, uint32_t, uint32_t, int);
int32_t BIndexCompareAt(BIndex*, BRGBinary*, int64_t, int64_t, int);
int32_t BIndexCompareRead(BIndex*, BRGBinary*, char*, int64_t, int);
uint32_t BIndexGetHashIndex(BIndex*, BRGBinary*, uint32_t, int);
uint32_t BIndexGetHashIndexFromRead(BIndex*, BRGBinary*, char*, int32_t, int);
void BIndexPrintReadMasked(BIndex*, char*, int, FILE*);
void BIndexInitialize(BIndex*);
#endif

