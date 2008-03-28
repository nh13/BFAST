#ifndef RGTREE_H_
#define RGTREE_H_

#include <stdio.h>
#include "BLibDefinitions.h"
#include "RGSeqPair.h"

int RGTreeInsert(RGTree*, char*, char*, int, int, int);
void RGTreeCleanUpTree(RGTree*);
void RGTreeInsertionSortNodes(RGTree*);
void RGTreeMergeSortNodes(RGTree*, int, int);
int RGTreeGetIndex(RGTree*, int, int);
void RGTreeDelete(RGTree*);
double RGTreeGetSize(RGTree*, int);
void RGTreePrintTree(FILE*, RGTree*);
int RGTreeReadTree(RGTree*, FILE*);
void RGTreePrintHeader(FILE*, RGTree*);
void RGTreeReadHeader(FILE*, RGTree*);
int RGTreeGetMatches(RGTree*, int, int, char, RGMatch*);
int GetIndexFromSequence(char*, int);
void RGNodeCopy(RGNode *src, RGNode *dest);
int RGNodeCompare(RGNode *a, RGNode *b);

#endif

