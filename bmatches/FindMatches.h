#ifndef FINDMATCHES_H_
#define FINDMATHCES_H_
#include <stdio.h>
#include "../blib/RGTree.h"

void RunMatches(char*, int, char*, char*, char*, char*, int, int, int, int, int, int, int, int, int, int);
int FindMatchesInIndexes(char**, int, int, int, FILE**, FILE*, int, int*, int*, int*);
void FindMatchesInIndex(FILE*, FILE*, RGIndex*, int);
int FindMatchesInTrees(char**, int, int, int*, int, int, int, int, int, int, int, FILE*, FILE*, int, int*, int*, int*);
void FindMatchesInTree(FILE*, FILE*, RGTree*, int**, int, int, int, int, int, int, int); 

#endif
