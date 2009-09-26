#ifndef RUNMATCH_H_ 
#define RUNMATCH_H_ 

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
#include "BLibDefinitions.h"

typedef struct {
	FILE *tempSeqFP;
	gzFile tempOutputFP;
	RGIndex *index;
	RGBinary *rg;
	int *offsets;
	int numOffsets;
	int space;
	int maxKeyMatches;
	int maxNumMatches;
	int whichStrand;
	int numMatches;
	int queueLength;
	int threadID;
} ThreadIndexData;

void FindMatches(
		char *fastaFileName,
		char *mainIndexes,
		char *secondaryIndexes,
		char *readFileName,
		char *offsets,
		int space,
		int startReadNum,
		int endReadNum,
		int keySize,
		int maxKeyMatches,
		int maxNumMatches,
		int whichStrand,
		int numThreads,
		int queueLength,
		char *tmpDir,
		int timing
		);
int FindMatchesInIndexes(char **indexFileNames,
		int32_t **indexIDs,
		int numRGIndexes,
		RGBinary *rg,
		int *offsets,
		int numOffsets,
		int colorSpace,
		int keySize,
		int maxKeyMatches,
		int maxNumMatches,
		int whichStrand,
		int numThreads,
		int queueLength,
		FILE ***tempSeqFPs,
		char ***tempSeqFileNames,
		gzFile outputFP,
		int copyForNextSearch,
		int indexesType,
		char *tmpDir,
		int timing,
		int *totalDataStructureTime,
		int *totalSearchTime,
		int *totalOutputTime);
int FindMatchesInIndex(char *indexFileName,
		RGBinary *rg,
		int *offsets,
		int numOffsets,
		int colorSpace,
		int keySize,
		int maxKeyMatches,
		int maxNumMatches,
		int whichStrand,
		int numThreads,
		int queueLength,
		FILE ***tempSeqFPs,
		gzFile indexFP,
		char *tmpDir,
		int timing,
		int *totalDataStructureTime,
		int *totalSearchTime,
		int *totalOutputTime);
void *FindMatchesInIndexThread(void *arg);

#endif
