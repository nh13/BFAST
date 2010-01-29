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
	gzFile tempRGMatchesFP;
	pthread_mutex_t *outputFP_mutex;
	gzFile outputFP;
	int32_t outputOffsets;
	int32_t *outputFP_threadID;
	int32_t numThreads;
	RGIndex *indexes;
	int32_t numIndexes;
	RGBinary *rg;
	int32_t *offsets;
	int numOffsets;
	int space;
	int maxKeyMatches;
	int maxNumMatches;
	int whichStrand;
	int numMatches;
	int queueLength;
	int threadID;
} ThreadIndexData;

void RunMatch(
		char *fastaFileName,
		char *mainIndexes,
		char *secondaryIndexes,
		char *readFileName,
		char *offsets,
		int loadAllIndexes,
		int compression,
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
		int timing,
		FILE *fpOut
		);
int FindMatchesInIndexSet(char **indexFileNames,
		int32_t **indexIDs,
		int numRGIndexes,
		RGBinary *rg,
		int32_t *offsets,
		int numOffsets,
		int loadAllIndexes,
		int colorSpace,
		int keySize,
		int maxKeyMatches,
		int maxNumMatches,
		int whichStrand,
		int numThreads,
		int queueLength,
		gzFile *tempRGMatchesFPs,
		char **tempRGMatchesFileNames,
		gzFile outputFP,
		int copyForNextSearch,
		int indexesType,
		char *tmpDir,
		int timing,
		int *totalDataStructureTime,
		int *totalSearchTime,
		int *totalOutputTime);
int FindMatches(char **indexFileName,
		int32_t numIndexes,
		RGBinary *rg,
		int32_t *offsets,
		int numOffsets,
		int loadAllIndexes,
		int colorSpace,
		int keySize,
		int maxKeyMatches,
		int maxNumMatches,
		int whichStrand,
		int numThreads,
		int queueLength,
		gzFile *tempRGMatchesFPs,
		char **tempRGMatchesFileNames,
		gzFile outputFP,
		int outputOffsets,
		char *tmpDir,
		int timing,
		int *totalDataStructureTime,
		int *totalSearchTime,
		int *totalOutputTime);
void *FindMatchesThread(void *arg);

#endif
