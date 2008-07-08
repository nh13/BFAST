#ifndef BINDEXDIST_H_
#define BINDEXDIST_H_

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

void PrintDistribution(RGIndex*, RGBinary*, char*, int, char*);
void GetMatchesFromChrPos(RGIndex*, RGBinary*, uint32_t, uint32_t, int, int64_t*, int64_t*, char*, char*);
void ReverseFile(FILE**, char**, char*);
void MergeFiles(FILE*, FILE*, FILE*);

#endif
