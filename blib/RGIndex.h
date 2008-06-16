#ifndef RGINDEX_H_
#define RGINDEX_H_

#include <stdio.h>
#include "RGBinary.h"
#include "BLibDefinitions.h"

void RGIndexCreate(RGIndex*, RGBinary*, RGIndexLayout*, int32_t, int32_t, int32_t, int32_t, int32_t, int32_t, int32_t, int32_t);
void RGIndexCreateHash(RGIndex*, RGBinary*);
void RGIndexSortNodes(RGIndex*, RGBinary*, int32_t);
void *RGIndexQuickSortNodes(void*);
void RGIndexQuickSortNodesHelper(RGIndex*, RGBinary*, int64_t, int64_t, int32_t);
void RGIndexQuickSortNodesGetPivots(RGIndex*, RGBinary*, int64_t, int64_t, int64_t*, int32_t, int32_t);
void RGIndexShellSortNodesHelper(RGIndex*, RGBinary*, int64_t, int64_t, int32_t, double*, int64_t, int64_t);
void RGIndexDelete(RGIndex*);
double RGIndexGetSize(RGIndex*, int32_t);
void RGIndexPrint(FILE*, RGIndex*, int32_t);
void RGIndexRead(FILE*, RGIndex*, int32_t);
void RGIndexPrintInfo(FILE*, int32_t);
void RGIndexPrintHeader(FILE*, RGIndex*, int32_t);
void RGIndexReadHeader(FILE*, RGIndex*, int32_t);
void RGIndexGetMatches(RGIndex*, RGBinary*, char*, int32_t, int8_t, int32_t, RGMatch*, int32_t);
int64_t RGIndexGetIndex(RGIndex*, RGBinary*, int64_t, int64_t, char*, int64_t*, int64_t*);
void RGIndexSwapAt(RGIndex*, int64_t, int64_t);
int64_t RGIndexGetPivot(RGIndex*, RGBinary*, int64_t, int64_t);
int32_t RGIndexCompareAt(RGIndex*, RGBinary*, int64_t, int64_t, int);
int32_t RGIndexCompareRead(RGIndex*, RGBinary*, char*, int64_t, int);
uint32_t RGIndexGetHashIndex(RGIndex*, RGBinary*, uint32_t, int);
uint32_t RGIndexGetHashIndexFromRead(RGIndex*, RGBinary*, char*, int32_t, int);
void RGIndexPrintReadMasked(RGIndex*, char*, FILE*);
#endif

