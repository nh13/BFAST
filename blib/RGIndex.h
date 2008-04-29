#ifndef RGINDEX_H_
#define RGINDEX_H_

#include <stdio.h>
#include "RGBinary.h"
#include "BLibDefinitions.h"

void RGIndexCreate(RGIndex*, RGBinary*, int, int);
void RGIndexCleanUpIndex(RGIndex*, RGBinary*, int);
void RGIndexSortNodes(RGIndex*, RGBinary*, int);
void *RGIndexQuickSortNodes(void*);
void RGIndexQuickSortNodesHelper(RGIndex*, RGBinary*, unsigned long long int, unsigned long long int, int, double*, unsigned long long int, unsigned long long int, int, unsigned long long int*, int, int);
void RGIndexDelete(RGIndex*);
double RGIndexGetSize(RGIndex*, int);
void RGIndexPrint(FILE*, RGIndex*, int);
int RGIndexRead(FILE*, RGIndex*, int);
void RGIndexPrintHeader(FILE*, RGIndex*, int);
void RGIndexReadHeader(FILE*, RGIndex*, int);
int RGIndexGetMatches(RGIndex*, RGBinary*, char*, char, int, RGMatch*);
long long int RGIndexGetFirstIndex(RGIndex*, RGBinary*, char*);
int RGIndexCompareAt(RGIndex*, RGBinary*, unsigned int, unsigned int);
int RGIndexCompareRead(RGIndex*, RGBinary*, char*, unsigned int);
#endif

