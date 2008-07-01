#ifndef BINDEXSTAT_H_
#define BINDEXSTAT_H_

#include "../blib/RGIndex.h"
#include "../blib/RGBinary.h"
#include "../blib/RGMatch.h"

typedef struct {
	int **counts;
	int *maxCount;
} Counts;

void PrintHistogram(RGIndex*, RGBinary*, int, int, char*);
void GetMatchesFromChrPos(RGIndex*, RGBinary*, uint32_t, uint32_t, int, int*, int*);

#endif
