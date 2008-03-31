#ifndef RGINDEX_H_
#define RGINDEX_H_

#include <stdio.h>
#include "BLibDefinitions.h"

int RGIndexInsert(RGIndex*, char*, int, int, int);
void RGIndexCleanUpIndex(RGIndex *index);
void RGIndexQuickSortNodes(RGIndex*, int, int, int); 
int RGIndexGetIndex(RGIndex*, unsigned char*);
void RGIndexDelete(RGIndex*);
double RGIndexGetSize(RGIndex*, int);
void RGIndexPrintIndex(FILE*, RGIndex*);
int RGIndexReadIndex(FILE*, RGIndex*);
void RGIndexPrintHeader(FILE*, RGIndex*);
void RGIndexReadHeader(FILE*, RGIndex*);
int RGIndexGetMatches(RGIndex*, unsigned char*, char, RGMatch*);
void RGIndexNodeCopy(RGIndexNode*, RGIndexNode*, int);
int RGIndexNodeCompare(RGIndexNode*, RGIndexNode*, int);
void RGIndexQuickSortNode(RGIndex*, int, int, int);
void RGIndexGetIndexFromSequence(char*, int, unsigned char*);
#endif

