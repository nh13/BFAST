#ifndef RGTREE_H_
#define RGTREE_H_

#include <stdio.h>
#include "BLibDefinitions.h"
#include "RGSeqPair.h"

int RGTreeInsert(RGTree*, char*, char*, int, int, int);
void RGTreeCleanUpTree(RGTree*);
void RGTreeQuickSortNodes(RGTree*, int, int, int);
int RGTreeGetIndex(RGTree*, int, int);
void RGTreeDelete(RGTree*);
double RGTreeGetSize(RGTree*, int);
void RGTreePrintTree(FILE*, RGTree*, int);
int RGTreeReadTree(FILE*, RGTree*, int);
void RGTreePrintHeader(FILE*, RGTree*, int);
void RGTreeReadHeader(FILE*, RGTree*, int);
int RGTreeGetMatches(RGTree*, int, int, char, int, RGMatch*);
int RGTreeGetIndexFromSequence(char*, int);
void RGTreeNodeCopy(RGTreeNode *src, RGTreeNode *dest);
int RGTreeNodeCompare(RGTreeNode *a, RGTreeNode *b);
void RGTreeQuickSortNode(RGTree*, int, int, int);

#endif
