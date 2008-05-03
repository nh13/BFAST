#ifndef FINDMATCHES_H_
#define FINDMATHCES_H_

/*
 *   _REENTRANT to grab thread-safe libraries
 *   _POSIX_SOURCE to get POSIX semantics
 */
#ifndef _REENTRANT
#define _REENTRANT
#endif
#ifndef _POSIX_SOURCE
#define _POSIX_SOURCE
#endif

#include <stdio.h>
#include "../blib/BLibDefinitions.h"

typedef struct {
	FILE *tempSeqFP;
	FILE *tempOutputFP;
	RGIndex *index;
	RGBinary *rg;
	int *offsets;
	int numOffsets;
	int numMismatches;
	int numInsertions;
	int numDeletions;
	int numGapInsertions;
	int numGapDeletions;
	int pairedEnd;
	int numMatches;
	int threadID;
} ThreadIndexData;

void FindMatches(char*, int, char*, char*, char*, char*, char*, int, int, int, int, int, int, int, int, int, int, int);
int FindMatchesInIndexes(char **rgIndexFileNames,
		int binaryInput,
		RGBinary *rg,
		int numRGIndexes,
		int *offsets,
		int numOffsets,
		int numMismatches,
		int numInsertions,
		int numDeletions,
		int numGapInsertions,
		int numGapDeletions,
		int pairedEnd,
		int numThreads,
		FILE ***tempSeqFPs,
		FILE *outputFP,
		int MainIndexes,
		int timing,
		int *totalDataStructureTime,
		int *totalSearchTime,
		int *totalOutputTime);
void *FindMatchesInIndex(void *arg);

#endif
