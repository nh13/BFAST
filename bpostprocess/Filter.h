#ifndef FILTER_H_
#define FILTER_H_

#include "../blib/AlignedRead.h"

int FilterAlignedRead(AlignedRead *a,
		int algorithm,
		int minScores,
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
		int unpaired,
		double mismatchScore,
		int avgMismatchQuality,
		int space);

int FilterAlignedEntry(AlignedEntry *a,
		int space,
		int minScores,
		int startChr,
		int startPos,
		int endChr,
		int endPos,
		int maxMismatches,
		int maxColorErrors);

void FilterReplaceAlignmentScore(AlignedEnd *a,
		double mismatchScore,
		int avgMismatchQuality,
		int space);

#endif
