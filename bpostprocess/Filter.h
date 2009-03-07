#ifndef FILTER_H_
#define FILTER_H_

#include "../blib/AlignedRead.h"

int FilterAlignedRead(AlignedRead *a,
		int algorithmReads,
		int mappingQuality,
		int minMappingQuality,
		int minScoreReads,
		int startChr,
		int startPos,
		int endChr,
		int endPos,
		int maxMismatches,
		int maxColorErrors,
		int minDistancePaired,
		int maxDistancePaired,
		int useDistancePaired,
		int contigAbPaired,
		int inversionsPaired,
		int unpaired);

int FilterAlignedEntry(AlignedEntry *a,
		int space,
		int minScoreReads,
		int startChr,
		int startPos,
		int endChr,
		int endPos,
		int maxMismatches,
		int maxColorErrors);

double GetMappingQuality(AlignedEnd *end,
		int32_t index);

#endif
