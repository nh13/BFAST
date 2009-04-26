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
	int binaryOutput;
	RGIndex *index;
	RGBinary *rg;
	int *offsets;
	int numOffsets;
	int space;
	int numMismatches;
	int numInsertions;
	int numDeletions;
	int numGapInsertions;
	int numGapDeletions;
	int maxKeyMatches;
	int maxNumMatches;
	int whichStrand;
	int numMatches;
	int threadID;
} ThreadIndexData;

void FindMatches(int, char*, char*, char*, char*, char*, int, int, int, int, int, int, int, int, int, int, int, int, int, int, char*, char*, char*, int);
int FindMatchesInIndexes(char **rgIndexFileNames,
		int binaryInput,
		RGBinary *rg,
		int numRGIndexes,
		int *offsets,
		int numOffsets,
		int colorSpace,
		int numMismatches,
		int numInsertions,
		int numDeletions,
		int numGapInsertions,
		int numGapDeletions,
		int keySize,
		int maxKeyMatches,
		int maxNumMatches,
		int whichStrand,
		int numThreads,
		FILE ***tempSeqFPs,
		char ***tempSeqFileNames,
		FILE *outputFP,
		int binaryOutput,
		int copyForNextSearch,
		int indexesType,
		char *tmpDir,
		int timing,
		int *totalDataStructureTime,
		int *totalSearchTime,
		int *totalOutputTime);
int FindMatchesInIndex(char *indexFileName,
		int binaryInput,
		RGBinary *rg,
		int *offsets,
		int numOffsets,
		int colorSpace,
		int numMismatches,
		int numInsertions,
		int numDeletions,
		int numGapInsertions,
		int numGapDeletions,
		int keySize,
		int maxKeyMatches,
		int maxNumMatches,
		int whichStrand,
		int numThreads,
		FILE ***tempSeqFPs,
		FILE *indexFP,
		int binaryOutput,
		char *tmpDir,
		int timing,
		int *totalDataStructureTime,
		int *totalSearchTime,
		int *totalOutputTime);
void *FindMatchesInIndexThread(void *arg);

#endif
