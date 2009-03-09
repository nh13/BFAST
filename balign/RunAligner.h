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

#include "../blib/RGIndexAccuracy.h"
#include "../blib/BLibDefinitions.h"
#include "Definitions.h"

typedef struct {
	FILE *inputFP;
	char *inputFileName;
	FILE *outputFP;
	char *outputFileName;
	FILE *notAlignedFP;
	char *notAlignedFileName;
	RGBinary *rg;
	int space;
	int offsetLength;
	int usePairedEndLength;
	int pairedEndLength;
	int mirroringType;
	int forceMirroring;
	int binaryInput;
	int binaryOutput;
	ScoringMatrix *sm;
	int alignmentType;
	int bestOnly;
	int64_t numLocalAlignments;
	QS qs;
	int threadID;
} ThreadData;

void RunAligner(RGBinary*, char*, char*, char*, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, char*, char*, char*, int, int*, int*);
void RunDynamicProgramming(FILE*, RGBinary*, char*, RGIndexAccuracyMismatchProfile*, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, char*, FILE*, FILE*, int, int*, int*);
void *RunDynamicProgrammingThread(void *);
#endif
