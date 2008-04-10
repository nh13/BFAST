#ifndef FINDMATCHES_H_
#define FINDMATHCES_H_

/*
 *   _REENTRANT to grab thread-safe libraries
 *   _POSIX_SOURCE to get POSIX semantics
 */
#define _REENTRANT
#define _POSIX_SOURCE

#include <stdio.h>
#include "../blib/RGTree.h"

typedef struct {
	FILE *tempSeqFP;
	FILE *tempOutputFP;
	RGIndex *index;
	int pairedEnd;
	int numMatches;
	int threadID;
} ThreadIndexData;

typedef struct {
	int temp;
	FILE *tempSeqFP;
	FILE *tempOutputFP;
	RGTree *tree;
	int *offsets;
	int numOffsets;
	int numMismatches;
	int numInsertions;
	int numDeletions;
	int numGapInsertions;
	int numGapDeletions;
	int pairedEnd;
	int threadID;
} ThreadTreeData;

void RunMatches(char*, int, char*, char*, char*, char*, int, int, int, int, int, int, int, int, int, int, int);
int FindMatchesInIndexes(char**, int, int, int, int, FILE***, FILE*, int, int*, int*, int*);
void *FindMatchesInIndex(void *arg);
int FindMatchesInTrees(char**, int, int, int*, int, int, int, int, int, int, int, int, FILE***, FILE*, int, int*, int*, int*);
void *FindMatchesInTree(void *arg);

#endif
