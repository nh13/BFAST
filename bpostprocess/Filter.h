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
		int pairedEnd,
		int algorithmReadsPaired,
		int minScoreReadsPaired,
		int minDistancePaired,
		int maxDistancePaired,
		int meanDistancePaired);

int FilterFirstReadInAlignEntries(AlignEntries *a,
		int algorithmReads,
		int minScoreReads,
		int startChr,
		int startPos,
		int endChr,
		int endPos);

int FilterAlignEntry(AlignEntry *a,
		int minScoreReads,
		int startChr,
		int startPos,
		int endChr,
		int endPos);

int FilterOneAlignEntries(AlignEntries *a,
		int startChr,
		int startPos,
		int endChr,
		int endPos,
		int algorithmReadsPaired,
		int minScoreReadsPaired,
		int minDistancePaired,
		int maxDistancePaired,
		int meanDistancePaired);

#endif
