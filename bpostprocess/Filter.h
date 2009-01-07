#ifndef FILTER_H_
#define FILTER_H_

#include "../blib/AlignEntries.h"

int FilterAlignEntries(AlignEntries *a,
		int algorithmReads,
		int minScoreReads,
		int startChr,
		int startPos,
		int endChr,
		int endPos,
		int maxMismatches,
		int maxColorErrors,
		int pairedEnd,
		int algorithmReadsPaired,
		int minScoreReadsPaired,
		int minDistancePaired,
		int maxDistancePaired,
		int maxMismatchesPaired,
		int maxColorErrorsPaired);

int FilterReadInAlignEntries(AlignEntries *a,
		int algorithmReads,
		int minScoreReads,
		int startChr,
		int startPos,
		int endChr,
		int endPos,
		int maxMismatches,
		int maxColorErrors,
		int which);

int FilterAlignEntry(AlignEntry *a,
		int space,
		int minScoreReads,
		int startChr,
		int startPos,
		int endChr,
		int endPos,
		int maxMismatches,
		int maxColorErrors,
		int *returnNumMismatches,
		int *returnNumColorErrors);

int FilterOneAlignEntries(AlignEntries *a,
		int startChr,
		int startPos,
		int endChr,
		int endPos,
		int minScoreReads,
		int maxColorErrors,
		int pairedEnd,
		int algorithmReadsPaired,
		int minScoreReadsPaired,
		int minDistancePaired,
		int maxDistancePaired,
		int maxMismatchesPaired,
		int maxColorErrorsPaired);

#endif
