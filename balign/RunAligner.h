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

#include "../blib/BLibDefinitions.h"
#include "../blib/AlignEntry.h"
#include "Definitions.h"

typedef struct {
	FILE *inputFP;
	char *inputFileName;
	FILE *outputFP;
	char *outputFileName;
	FILE *notAlignedFP;
	char *notAlignedFileName;
	RGBinary *rgBinary;
	int colorSpace;
	int offsetLength;
	int pairedEnd;
	int usePairedEndLength;
	int pairedEndLength;
	int forceMirroring;
	int binaryInput;
	ScoringMatrix *sm;
	int threadID;
} ThreadData;

void RunAligner(RGBinary*, char*, char*, int, int, int, int, int, int, int, int, int, int, int, int, int, char*, char*, char*, int*, int*);
void RunDynamicProgramming(FILE*, RGBinary*, char*, int, int, int, int, int, int, int, int, int, int, int, int, int, char*, FILE*, FILE*, int*, int*);
void *RunDynamicProgrammingThread(void *);
void RunDynamicProgrammingThreadHelper(RGBinary*, uint8_t, uint32_t, int8_t, char*, int, int, int, ScoringMatrix*, AlignEntry*);
void GetSequenceFromReferenceGenome(RGBinary*, int, int, char, int, char*, int, int*, int*);
void GetReverseComplimentAnyCase(char*, char*, int);
#endif
