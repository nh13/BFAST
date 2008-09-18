#ifndef BINDEXHIST_H_
#define BINDEXHIST_H_

#include "../blib/RGIndex.h"
#include "../blib/RGBinary.h"
#include "../blib/RGMatch.h"

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
	int numDifferent;
	int64_t totalForward;
	int64_t totalReverse;
	int threadID;
} ThreadData;

void PrintHistogram(RGIndex*, RGBinary*, int, int, int, char*);
void *PrintHistogramThread(void *arg);
void GetPivots(RGIndex*, RGBinary*, int64_t*, int64_t*, int64_t);
void GetMatchesFromChrPos(RGIndex*, RGBinary*, uint32_t, uint32_t, int, int64_t*, int64_t*);

#endif
