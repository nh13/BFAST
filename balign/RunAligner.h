#ifndef RUNALIGNER_H_
#define RUNALIGNER_H_

/*
 *     _REENTRANT to grab thread-safe libraries
 *      _POSIX_SOURCE to get POSIX semantics
 */
#define _REENTRANT
#define _POSIX_SOURCE

#include "../blib/BLibDefinitions.h"
#include "Definitions.h"

typedef struct {
	FILE *inputFP;
	FILE *outputFP;
	RGBinary *rgBinary;
	int offsetLength;
	int maxNumMatches;
	int pairedEnd;
	ScoringMatrix *sm;
	int threadID;
} ThreadData;

void RunAligner(RGBinary*, char*, char*, int, int, int, int, int, char*, char*);
void RunDynamicProgramming(FILE*, RGBinary*, char*, int, int, int, int, FILE*);
void *RunDynamicProgrammingThread(void *);
void GetSequenceFromReferenceGenome(RGBinary*, int, int, char, int, char*, int, int*, int*);
void GetReverseComplimentAnyCase(char*, char*, int);
#endif
