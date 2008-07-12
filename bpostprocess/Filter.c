#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <limits.h>
#include "../blib/BError.h"
#include "../blib/AlignEntries.h"
#include "Definitions.h"
#include "Filter.h"

/* TODO */
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
		int meanDistancePaired)
{
	int foundType=NoneFound;
	assert(a->pairedEnd == pairedEnd);

	/* We should only modify "a" if it is going to be reported */ 

	/* Filter for reads that are not paired end */
	if(pairedEnd == 0) {
		foundType=FilterFirstReadInAlignEntries(a,
				algorithmReads,
				minScoreReads,
				startChr,
				startPos,
				endChr,
				endPos);
	}
	else {
		/* We should return whether we printed as an chrAb or inversion */
		foundType=FilterOneAlignEntries(a,
				startChr,
				startPos,
				endChr,
				endPos,
				algorithmReadsPaired,
				minScoreReadsPaired,
				minDistancePaired,
				maxDistancePaired,
				meanDistancePaired);
	}

	return foundType;
}

/* TODO */
int FilterFirstReadInAlignEntries(AlignEntries *a,
		int algorithmReads,
		int minScoreReads,
		int startChr,
		int startPos,
		int endChr,
		int endPos)
{
	int i;
	int numNotFiltered = 0;
	int uniqueIndex=-1;
	int bestScore = INT_MAX;
	double numBestScore = 0;
	int bestScoreIndex = -1;
	int foundType = Found;

	/* Filter each entry for this read */
	for(i=0;i<a->numEntriesOne && Found==foundType;i++) {
		/* Check if we should filter */
		if(1!=FilterAlignEntry(&a->entriesOne[i],
					minScoreReads,
					startChr,
					startPos,
					endChr,
					endPos)) {
			uniqueIndex=i;
			numNotFiltered++;
			if(numNotFiltered > 1) {
				/* This mean that two ore more entries passed the filters.  Therefore
				 * we will not have a unique alignment for this read.  Free memory and 
				 * return */
				foundType = NoneFound;
			}
			else if(1==algorithmReads) {
				if(a->entriesOne[i].score == bestScore) {
					assert(numBestScore > 0);
					numBestScore++;
				}
				else if(a->entriesOne[i].score > bestScore) {
					bestScore = a->entriesOne[i].score;
					numBestScore = 1;
					bestScoreIndex = i;
				}
			}
		}
	}

	/* Check if we found an alignment */
	if(0==algorithmReads && Found==foundType) {
		assert(uniqueIndex >= 0 && uniqueIndex < a->numEntriesOne);
		/* Copy to the front and update */
		AlignEntriesKeepOnly(a, uniqueIndex, 0, 0);
		foundType=Found;
	}
	else if(1==algorithmReads && numBestScore == 1) {
		assert(foundType == Found);
		/* Check if there is a best score (unique best score) */
		AlignEntriesKeepOnly(a, bestScoreIndex, 0, 0);
		foundType = Found;
	}
	else {
		foundType = NoneFound;
	}

	return foundType;
}

/* TODO */
int FilterAlignEntry(AlignEntry *a,
		int minScoreReads,
		int startChr,
		int startPos,
		int endChr,
		int endPos)
{
	/* Check if the alignment has at least the minimum score */
	if(a->score < minScoreReads) {
		return 1;
	}
	/* Check the genomic location */
	if(a->chromosome < startChr || 
			(a->chromosome == startChr && a->position < startPos) ||
			(a->chromosome == endChr && a->position > endPos) ||
			(a->chromosome > endChr)) {
		return 1;
	}
	return 0;
}

/* TODO */
int FilterOneAlignEntries(AlignEntries *a,
		int startChr,
		int startPos,
		int endChr,
		int endPos,
		int algorithmReadsPaired,
		int minScoreReadsPaired,
		int minDistancePaired,
		int maxDistancePaired,
		int meanDistancePaired)
{
	char *FnName = "FilterOneAlignEntries";
	/* best score */
	int bestScoreInitialized=0;
	int bestScoreOne=-1, bestScoreTwo=-1;
	int numBestScore=0;
	/* unique */
	int uniqueOne=-1, uniqueTwo=-1;
	int numUnique=0;
	/* closest to the mean unique/best score */
	int closestToMeanOne=-1, closestToMeanTwo=-1;
	int closestToMean=INT_MAX;
	int curDistanceToMean = -1;
	int closestToMeanScore=INT_MIN;
	int numClosestToMean=0;
	/* best score outside */
	int bestScoreOutsideInitialized=0;
	int bestScoreOutsideOne=-1, bestScoreOutsideTwo=-1;
	int numBestScoreOutside=0;
	/* unique outside */
	int uniqueOutsideOne=-1, uniqueOutsideTwo=-1;
	int numUniqueOutside=0;
	/* current distance to the mean */
	int curChrDistance, curPositionDistance;
	int i, j;
	int foundType=NoneFound;

	/* Note: this does not take into account strandedness */

	/* Go through all possible pairs of alignments */
	for(i=0;i<a->numEntriesOne;i++) { /* First read */
		for(j=0;j<a->numEntriesTwo;j++) { /* Second read */
			/* Get the distance between the two */
			curChrDistance = a->entriesTwo[j].chromosome - a->entriesOne[i].chromosome;
			curPositionDistance = a->entriesTwo[j].position - a->entriesOne[i].position;
			if(minScoreReadsPaired <= (a->entriesOne[i].score + a->entriesTwo[j].score)) {
				if(0==curChrDistance &&
						curPositionDistance <= maxDistancePaired &&
						curPositionDistance >= minDistancePaired) {
					/* Inside the bounds */

					/* best score */
					if(0==bestScoreInitialized ||
							(a->entriesOne[bestScoreOne].score + a->entriesTwo[bestScoreTwo].score) <=
							(a->entriesOne[i].score + a->entriesTwo[j].score)) {
						if((a->entriesOne[bestScoreOne].score + a->entriesTwo[bestScoreTwo].score) ==
								(a->entriesOne[i].score + a->entriesTwo[j].score)) {
							/* If we have the same score update the number */ 
							assert(numBestScore > 0);
							numBestScore++;
						}
						else {
							/* Update the score */
							bestScoreOne = i;
							bestScoreTwo = j;
							numBestScore=1;
							bestScoreInitialized = 1;
						}
					}
					/* unique */
					numUnique++;
					uniqueOne = i;
					uniqueTwo = j;
					/* closest to the mean unique/best score */ 
					curDistanceToMean = abs(curPositionDistance - meanDistancePaired);
					if(algorithmReadsPaired == 2) { /* closest unique */
						if(curDistanceToMean == closestToMean) {
							/* Increment */
							numClosestToMean++;
						}
						else if(curDistanceToMean < closestToMean) {
							/* Update */
							closestToMean = curDistanceToMean;
							closestToMeanOne = i;
							closestToMeanTwo = j;
							numClosestToMean = 1;
						}
					}
					else if(algorithmReadsPaired == 3) { /* closest best score */
						if(curDistanceToMean == closestToMean) {
							if((a->entriesOne[closestToMeanOne].score + a->entriesTwo[closestToMeanTwo].score) < 
									(a->entriesOne[i].score + a->entriesTwo[j].score)) {
								/* Update */
								closestToMean = curDistanceToMean;
								closestToMeanOne = i;
								closestToMeanTwo = j;
								closestToMeanScore = (a->entriesOne[i].score + a->entriesTwo[j].score);
								numClosestToMean = 1;
							}
							else if((a->entriesOne[closestToMeanOne].score + a->entriesTwo[closestToMeanTwo].score) == 
									(a->entriesOne[i].score + a->entriesTwo[j].score)) {
								numClosestToMean++;
							}
						}
						else if(curDistanceToMean < closestToMean) {
							/* Update */
							closestToMean = curDistanceToMean;
							closestToMeanOne = i;
							closestToMeanTwo = j;
							closestToMeanScore = (a->entriesOne[i].score + a->entriesTwo[j].score);
							numClosestToMean = 1;
						}
					}
				}
				else {
					/* Outside the bounds */

					/* best score outside */
					if(0==bestScoreOutsideInitialized ||
							(a->entriesOne[bestScoreOutsideOne].score + a->entriesTwo[bestScoreOutsideTwo].score) <=
							(a->entriesOne[i].score + a->entriesTwo[j].score)) {
						if((a->entriesOne[bestScoreOutsideOne].score + a->entriesTwo[bestScoreOutsideTwo].score) ==
								(a->entriesOne[i].score + a->entriesTwo[j].score)) {
							/* If we have the same score update the number */ 
							assert(numBestScoreOutside > 0);
							numBestScoreOutside++;
						}
						else {
							/* Update the score */
							bestScoreOutsideOne = i;
							bestScoreOutsideTwo = j;
							numBestScoreOutside=1;
							bestScoreOutsideInitialized = 1;
						}
					}
					/* unique outside */
					numUniqueOutside++;
					uniqueOutsideOne = i;
					uniqueOutsideTwo = j;
				}
			}
		}
	}

	/* Check if we have a unique pair of alignments */
	switch(algorithmReadsPaired) {
		case 0: /* unique */
			if(numUnique == 1) { 
				/* Prefer within the boundaries */
				AlignEntriesKeepOnly(a,
						uniqueOne,
						uniqueTwo,
						1);
				foundType=Found;
			}
			else if(numUnique == 0 &&
					numUniqueOutside == 1) { 
				/* If there are no pairs within the bounds */
				AlignEntriesKeepOnly(a,
						uniqueOutsideOne,
						uniqueOutsideTwo,
						1);
				foundType=OutsideBounds;
			}
			else {
				foundType=NoneFound;
			}
			break;
		case 1: /* best score */
			if(numBestScore == 1) {
				/* Prefer within the boundaries */
				AlignEntriesKeepOnly(a,
						bestScoreOne,
						bestScoreTwo,
						1);
				foundType=Found;
			}
			else if(numBestScore == 0 &&
					numBestScoreOutside == 1) {
				/* If there are no pairs within the bounds */
				AlignEntriesKeepOnly(a,
						bestScoreOutsideOne,
						bestScoreOutsideTwo,
						1);
				foundType=OutsideBounds;
			}
			else {
				foundType=NoneFound;
			}
			break;
		case 2: /* closest to the mean unique */
		case 3: /* closest to the mean best score */
			if(numClosestToMean == 1) {
				AlignEntriesKeepOnly(a,
						closestToMeanOne,
						closestToMeanTwo,
						1);
				foundType=Found;
			}
			else {
				foundType=NoneFound;
			}
			break;
		default:
			PrintError(FnName,
					"algorithmReadsPaired",
					"Could not understand algorithmReadsPaired",
					Exit,
					OutOfRange);
			break;
	}

	/* If outside bounds, print accordingly */
	if(OutsideBounds == foundType) {
		if(a->entriesOne[0].strand == a->entriesTwo[0].strand) {
			/* Same strand - chromosomal abnormality */
			foundType = ChrAb;
		}
		else {
			/* Different strand - inversion */
			foundType = Inversion;
		}
	}
	return foundType;
}
