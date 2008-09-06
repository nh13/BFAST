#ifndef BINDEXDIST_H_
#define BINDEXDIST_H_

#include "../blib/RGIndex.h"
#include "../blib/RGBinary.h"
#include "../blib/RGMatch.h"

typedef struct {
	char **reads;
	int64_t *readCounts;
	int64_t low;
	int64_t high;
	char *tmpDir;
	int readLength;
	int showPercentComplete;
	int threadID;
	int32_t colorSpace;
} ThreadData;

void PrintDistribution(RGIndex*, RGBinary*, char*, int, char*, int, int32_t);
void GetMatchesFromChrPos(RGIndex*, RGBinary*, uint32_t, uint32_t, int, int64_t*, int64_t*, char*, char*, int32_t);
void *MergeSortReads(void *arg);
void MergeSortReadsHelper(char**, int64_t*, int64_t, int64_t, int64_t, int64_t, int, double*, char*, int);
void MergeHelper(char**, int64_t*, int64_t, int64_t, int64_t, char*, int);

#endif
