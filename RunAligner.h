#ifndef RUNALIGNER_H_
#define RUNALIGNER_H_

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
	gzFile inputFP;
	char *inputFileName;
	gzFile outputFP;
	char *outputFileName;
	gzFile notAlignedFP;
	char *notAlignedFileName;
	RGBinary *rg;
	int space;
	int offsetLength;
	int usePairedEndLength;
	int pairedEndLength;
	int mirroringType;
	int forceMirroring;
	ScoringMatrix *sm;
	int alignmentType;
	int bestOnly;
	int64_t numLocalAlignments;
	int32_t avgMismatchQuality;
	double mismatchScore;
	int queueLength;
	int threadID;
} ThreadData;

void RunAligner(RGBinary*, char*, char*, int, int, int, int, int, int, int, int, int, int, int, int, int, int, char*, char*, char*, int*, int*);
void RunDynamicProgramming(gzFile, RGBinary*, char*, int, int, int, int, int, int, int, int, int, int, int, int, int, int, char*, gzFile, gzFile, int*, int*);
void *RunDynamicProgrammingThread(void *);
int32_t GetMatches(gzFile, RGMatches*, int32_t);
#endif
