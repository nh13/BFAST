#ifndef RGTREE_H_
#define RGTREE_H_

#include <stdio.h>
#include "BLibDefinitions.h"
#include "RGSeqPair.h"

int RGTreeInsert(RGTree*, char*, char*, int, int, int);
void RGTreeDelete(RGTree*);
void RGTreeDeleteHelper(void**, int, int);
double RGTreeGetSize(RGTree*, int);
double RGTreeGetSizeHelper(void*, int, int, int);
int GetIndexFromSequence(char*, int);
int RGTreeGetMatches(RGTree*, int, int, char, RGMatch*);

void RGTreePrintTree(FILE*, RGTree*);
void RGTreePrintTreeHelper(FILE*, void*, int, int);
int RGTreeReadFromFile(RGTree*, FILE*);
int RGTreeReadFromFileHelper(RGNode*, FILE*);
void RGTreePrintHeader(FILE*, RGTree*);
void RGTreeReadHeader(FILE*, RGTree*);

#endif

