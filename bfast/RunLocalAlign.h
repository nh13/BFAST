#ifndef RUNLOCALALIGN_H_
#define RUNLOCALALIGN_H_

/*
 *     _REENTRANT to grab thread-safe libraries
 *      _POSIX_SOURCE to get POSIX semantics
 */
#ifndef _REENTRANT
#define _REENTRANT
#endif
#ifndef _POSIX_SOURCE
#define _POSIX_SOURCE
#endif

#include "BLibDefinitions.h"

typedef struct {
	pthread_mutex_t matchFP_mutex;
	gzFile matchFP;
	int32_t matchFPctr;
	gzFile outputFP;
	pthread_mutex_t outputFP_mutex;
	int32_t outputFPctr;
	// No mutex for these
	int32_t startReadNum;
	int32_t endReadNum;
} ThreadFileData;

typedef struct {
	ThreadFileData *fdata;
	RGBinary *rg;
	int32_t space;
	int32_t offsetLength;
	int32_t usePairedEndLength;
	int32_t pairedEndLength;
	int32_t mirroringType;
	int32_t forceMirroring;
	ScoringMatrix *sm;
	int32_t ungapped;
	int32_t unconstrained;
	int32_t bestOnly;
	int64_t numLocalAlignments;
	int32_t avgMismatchQuality;
	double mismatchScore;
	int32_t queueLength;
	int32_t threadID;
	int64_t numAligned;
	int64_t numNotAligned;
	int32_t fileTime;
} ThreadData;

void RunAligner(char*, char*, char*, int32_t, int32_t, int32_t, int32_t, int32_t, int32_t, int32_t, int32_t, int32_t, int32_t, int32_t, int32_t, int32_t, int32_t, int32_t, int32_t*, int32_t*, int32_t*);
void RunDynamicProgramming(gzFile, RGBinary*, char*, int32_t, int32_t, int32_t, int32_t, int32_t, int32_t, int32_t, int32_t, int32_t, int32_t, int32_t, int32_t, int32_t, int32_t, int32_t, gzFile, int32_t*, int32_t*);
void *RunDynamicProgrammingThread(void *);
int32_t GetMatches(ThreadFileData*, RGMatches*, int32_t, int32_t*);
void SkipMatches(ThreadFileData*);
#endif
