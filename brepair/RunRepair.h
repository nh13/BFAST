#ifndef RUNREPAIR_H_
#define RUNREPAIR_H_

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

#include "../blib/BLibDefinitions.h"

typedef struct {
	FILE *inputFP;
	char *inputFileName;
	FILE *outputFP;
	char *outputFileName;
	FILE *notRepairedFP;
	char *notRepairedFileName;
	RGBinary *rg;
	int space;
	int minPairedEndDistance;
	int maxPairedEndDistance;
	int mirroringType;
	int strandedness;
	int binaryInput;
	int binaryOutput;
	ScoringMatrix *sm;
	int alignmentType;
	int64_t numLocalAlignments;
	int threadID;
} RepairThreadData;

void RunRepair(RGBinary*, char*, char*, int, int, int, int, int, int, int, int, char*, char*, char*, int, int*, int*);
void RunRepairHelper(FILE*, RGBinary*, char*, int, int, int, int, int, int, int, int, char*, FILE*, FILE*, int, int*, int*);
void *RunRepairHelperThread(void*);
#endif
