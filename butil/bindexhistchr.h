#ifndef BINDEXHISTCHR_H_
#define BINDEXHISTCHR_H_

#include "../bfast/RGIndex.h"
#include "../bfast/RGBinary.h"
#include "../bfast/RGMatch.h"

typedef struct {
	int64_t **counts;
	int64_t *maxCount;
} Counts;

typedef struct {
	int64_t startIndex;
	int64_t endIndex;
	RGIndex *index;
	RGBinary *rg;
	Counts c;
	int numMismatchesStart;
	int numMismatchesEnd;
	int whichStrand;
	int chr;
	int numDifferent;
	int64_t totalForward;
	int64_t totalReverse;
	int threadID;
} ThreadData;

void PrintHistogram(RGIndex*, RGBinary*, int, int, int, int, int, char*);
void *PrintHistogramThread(void *arg);
void GetPivots(RGIndex*, RGBinary*, int64_t*, int64_t*, int64_t);
int GetMatchesFromContigPos(RGIndex*, RGBinary*, uint32_t, uint32_t, int, int64_t*, int64_t*);

#endif
