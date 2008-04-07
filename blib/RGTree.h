#ifndef RGTREE_H_
#define RGTREE_H_

#include <stdio.h>
#include "BLibDefinitions.h"
#include "RGSeqPair.h"

int RGTreeInsert(RGTree*, char*, char*, unsigned int, unsigned int, unsigned int);
void RGTreeCleanUpTree(RGTree*);
void RGTreeQuickSortNodes(RGTree*, unsigned int, unsigned int, unsigned int);
unsigned int RGTreeGetIndex(RGTree*, unsigned int, unsigned int);
void RGTreeDelete(RGTree*);
double RGTreeGetSize(RGTree*, int);
void RGTreePrintTree(FILE*, RGTree*, int);
int RGTreeReadTree(FILE*, RGTree*, int);
void RGTreePrintHeader(FILE*, RGTree*, int);
void RGTreeReadHeader(FILE*, RGTree*, int);
int RGTreeGetMatches(RGTree*, unsigned int, unsigned int, char, int, RGMatch*);
unsigned int RGTreeGetIndexFromSequence(char*, int);
void RGTreeNodeCopy(RGTreeNode *src, RGTreeNode *dest);
int RGTreeNodeCompare(RGTreeNode *a, RGTreeNode *b);
void RGTreeQuickSortNode(RGTree*, unsigned int, unsigned int, unsigned int);

#endif
