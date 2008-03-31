#ifndef FINDMATCHES_H_
#define FINDMATHCES_H_
#include <stdio.h>
#include "../blib/RGTree.h"

void RunMatches(char*, char*, char*, char*, char*, int, int, int, int, int, int);
int FindMatchesInIndexes(char**, int, int, FILE**, FILE*);
void FindMatchesInIndex(FILE*, FILE*, RGIndex*, int);
int FindMatchesInTrees(char**, int, int*, int, int, int, int, int, FILE*, FILE*);
void FindMatchesInTree(FILE*, FILE*, RGTree*, int**, int, int, int, int, int); 

#endif
