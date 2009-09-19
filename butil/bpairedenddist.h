#ifndef BPAIREDENDDIST_H_
#define BPAIREDENDDIST_H_

#include <zlib.h>
#include "../bfast/RGIndex.h"
#include "../bfast/RGBinary.h"
#include "../bfast/RGMatch.h"

typedef struct {
	int minDistance;
	int maxDistance;
	int binSize;
	int32_t *counts;
	int32_t numCounts;
} Bins;

void PrintDistributionFromBMF(gzFile, Bins*);
void PrintDistributionFromBAF(gzFile, Bins*);

int BinsInsert(Bins*, int32_t);
void BinsPrint(Bins*, FILE*);
#endif
