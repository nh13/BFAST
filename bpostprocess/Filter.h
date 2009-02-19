#ifndef FILTER_H_
#define FILTER_H_

#include "../blib/AlignEntries.h"

int FilterAlignEntries(AlignEntries *a,
		int algorithmReads,
		int uniquenessScore,
		int minUniquenessScore,
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
		int useDistancePaired,
		int contigAbPaired,
		int inversionsPaired,
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

double GetUniquenessScore(AlignEntry *a,
		int32_t numEntries,
		int32_t bestScoreIndex);

#endif
