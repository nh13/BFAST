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
		int minDistancePaired,
		int maxDistancePaired,
		int unpaired);

int FilterAlignEntry(AlignEntry *a,
		int space,
		int minScoreReads,
		int startChr,
		int startPos,
		int endChr,
		int endPos,
		int maxMismatches,
		int maxColorErrors);

int GetNumMismatchesInAlignEntry(AlignEntry *a);

int GetNumColorErrorsInAlignEntry(AlignEntry *a, int space);

#endif
