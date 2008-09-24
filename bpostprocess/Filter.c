#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <limits.h>
#include "../blib/BError.h"
#include "../blib/AlignEntry.h"
#include "../blib/AlignEntries.h"
#include "Definitions.h"
#include "Filter.h"

/* TODO */
int FilterAlignEntries(AlignEntries *a,
		int algorithmReads,
		int minScoreReads,
		int startContig,
		int startPos,
		int endContig,
		int endPos,
		int pairedEnd,
		int algorithmReadsPaired,
		int minScoreReadsPaired,
		int minDistancePaired,
		int maxDistancePaired,
		int meanDistancePaired)
{
	char *FnName="FilterAlignEntries";
	int foundType=NoneFound;
	AlignEntries tmpA;
	assert(a->pairedEnd == pairedEnd);
	AlignEntriesInitialize(&tmpA);

	/* We should only modify "a" if it is going to be reported */ 

	/* Filter for reads that are not paired end */
	if(pairedEnd == 0) {
		switch(algorithmReads) {
			case NoFiltering:
				/* Do nothing */
				foundType=Found;
				break;
			case AllNotFiltered:
				AlignEntriesCopy(a, &tmpA);
				foundType=FilterReadInAlignEntries(a,
						algorithmReads,
						minScoreReads,
						startContig,
						startPos,
						endContig,
						endPos,
						First);
				/* If we filtered all, copy back */
				if(foundType == NoneFound) {
					AlignEntriesCopy(&tmpA, a);
				}
				break;
			case Unique:
			case BestScore:
				foundType=FilterReadInAlignEntries(a,
						algorithmReads,
						minScoreReads,
						startContig,
						startPos,
						endContig,
						endPos,
						First);
				break;
			default:
				PrintError(FnName,
						"algorithmReads",
						"Could not understand algorithmReads",
						Exit,
						OutOfRange);
				break;
		}
	}
	else {
		switch(algorithmReadsPaired) {
			case NoFiltering:
				/* Do nothing */
				foundType=Found;
				break;
			case AllNotFiltered:
				/* Filter both ends separately */
				AlignEntriesCopy(a, &tmpA);
				if(Found==FilterReadInAlignEntries(a,
							algorithmReads,
							minScoreReads,
							startContig,
							startPos,
							endContig,
							endPos,
							First) &&
						Found==FilterReadInAlignEntries(a,
							algorithmReads,
							minScoreReads,
							startContig,
							startPos,
							endContig,
							endPos,
							Second)) {
					foundType=Found;
				}
				else {
				/* If we filtered all, copy back */
					AlignEntriesCopy(&tmpA, a);
				}
			case Unique:
			case BestScore:
			case MeanUnique:
			case MeanBestScore:
				foundType=FilterOneAlignEntries(a,
						startContig,
						startPos,
						endContig,
						endPos,
						algorithmReadsPaired,
						minScoreReadsPaired,
						minDistancePaired,
						maxDistancePaired,
						meanDistancePaired);
				break;
			default:
				PrintError(FnName,
						"algorithmReadsPaired",
						"Could not understand algorithmReads",
						Exit,
						OutOfRange);
				break;
		}
	}

	return foundType;
}

/* TODO */
int FilterReadInAlignEntries(AlignEntries *a,
		int algorithmReads,
		int minScoreReads,
		int startContig,
		int startPos,
		int endContig,
		int endPos,
		int which)
{
	char *FnName="FilterReadInAlignEntries";
	int i;
	int numNotFiltered = 0;
	int uniqueIndex=-1;
	int bestScore = INT_MIN;
	double numBestScore = 0;
	int bestScoreIndex = -1;
	int foundType = Found;
	AlignEntry **ptrEntries=NULL;
	int *ptrNumEntries=NULL;

	switch(which) {
		case First:
			ptrEntries = &a->entriesOne;
			ptrNumEntries = &a->numEntriesOne;
			break;
		case Second:
			ptrEntries = &a->entriesTwo;
			ptrNumEntries = &a->numEntriesTwo;
			break;
		default:
			PrintError(FnName,
					"which",
					"Could not understand which",
					Exit,
					OutOfRange);
			break;
	}

	/* Filter each entry for this read */
	for(i=0;i<(*ptrNumEntries) && Found==foundType;i++) {
		/* Check if we should filter */
		if(1!=FilterAlignEntry(&(*ptrEntries)[i],
					minScoreReads,
					startContig,
					startPos,
					endContig,
					endPos)) {
			uniqueIndex=i;
			numNotFiltered++;

			switch(algorithmReads) {
				case AllNotFiltered:
					/* Do nothing */
					break;
				case Unique:
					if(numNotFiltered > 1) {
						/* This mean that two ore more entries passed the filters.  Therefore
						 * we will not have a unique alignment for this read.  Free memory and 
						 * return */
						foundType = NoneFound;
					}
					break;
				case BestScore:
					if((*ptrEntries)[i].score == bestScore) {
						assert(numBestScore > 0);
						numBestScore++;
					}
					else if((*ptrEntries)[i].score > bestScore) {
						bestScore = (*ptrEntries)[i].score;
						numBestScore = 1;
						bestScoreIndex = i;
					}
					break;
				default:
					PrintError(FnName,
							"algorithmReads",
							"Could not understand algorithmReads",
							Exit,
							OutOfRange);
					break;
			}
		}
		else {
			switch(algorithmReads) {
				case AllNotFiltered:
					/* Copy end here */
					if(i < (*ptrNumEntries)-1) {
						AlignEntryCopy(&(*ptrEntries)[(*ptrNumEntries)-1], &(*ptrEntries)[i]);
					}
					/* Reallocate */
					switch(which) {
						case First:
							AlignEntriesReallocate(a, a->numEntriesOne-1, a->numEntriesTwo, a->pairedEnd, a->space); 
							break;
						case Second:
							AlignEntriesReallocate(a, a->numEntriesOne, a->numEntriesTwo-1, a->pairedEnd, a->space); 
							break;
						default:
							PrintError(FnName,
									"which",
									"Could not understand which",
									Exit,
									OutOfRange);
							break;
					}
					/* Since we removed, do not incrememt */
					i--;
					break;
				case Unique:
				case BestScore:
					/* Do nothing */
					break;
				default:
					PrintError(FnName,
							"algorithmReads",
							"Could not understand algorithmReads",
							Exit,
							OutOfRange);
					break;
			}
		}
	}

	/* Check if we found an alignment */
	switch(algorithmReads) {
		case AllNotFiltered:
			if((*ptrNumEntries) > 0) {
				foundType = Found;
			}
			else {
				foundType = NoneFound;
			}
			break;
		case Unique:
			/* Do nothing */
			if(foundType == Found) {
				/* Copy to the front and update */
				switch(which) {
					case First:
						AlignEntriesKeepOnly(a, uniqueIndex, 0, 0, a->space);
						break;
					case Second:
						AlignEntriesKeepOnly(a, 0, uniqueIndex, 0, a->space);
						break;
					default:
						PrintError(FnName,
								"which",
								"Could not understand which",
								Exit,
								OutOfRange);
						break;
				}
			}
			else {
				foundType = NoneFound;
			}
			break;
		case BestScore:
			if(numBestScore == 1) {
				assert(foundType == Found);
				/* Check if there is a best score (unique best score) */
				switch(which) {
					case First:
						AlignEntriesKeepOnly(a, bestScoreIndex, 0, 0, a->space);
						break;
					case Second:
						AlignEntriesKeepOnly(a, 0, bestScoreIndex, 0, a->space);
						break;
					default:
						PrintError(FnName,
								"which",
								"Could not understand which",
								Exit,
								OutOfRange);
						break;
				}
			}
			else {
				foundType = NoneFound;
			}
			break;
		default:
			PrintError(FnName,
					"algorithmReads",
					"Could not understand algorithmReads",
					Exit,
					OutOfRange);
			break;
	}
	return foundType;
}

/* TODO */
int FilterAlignEntry(AlignEntry *a,
		int minScoreReads,
		int startContig,
		int startPos,
		int endContig,
		int endPos)
{
	/* Check if the alignment has at least the minimum score */
	if(a->score < minScoreReads) {
		return 1;
	}
	/* Check the genomic location */
	if((startContig != 0 && a->contig < startContig) || 
			(startContig != 0 && a->contig == startContig && startPos != 0 && a->position < startPos) ||
			(endContig != 0 && a->contig == endContig && endPos != 0 && a->position > endPos) ||
			(endContig != 0 && a->contig > endContig)) {
		return 1;
	}
	return 0;
}

/* TODO */
int FilterOneAlignEntries(AlignEntries *a,
		int startContig,
		int startPos,
		int endContig,
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
	int curContigDistance, curPositionDistance;
	int i, j;
	int foundType=NoneFound;

	/* Note: this does not take into account strandedness */

	/* Go through all possible pairs of alignments */
	for(i=0;i<a->numEntriesOne;i++) { /* First read */
		/* Filter first entry for contig position */
		if(0 == FilterAlignEntry(&a->entriesOne[i],
					INT_MIN,
					startContig,
					startPos,
					endContig,
					endPos)) {
			for(j=0;j<a->numEntriesTwo;j++) { /* Second read */
				/* Filter second entry for contig position.  Also filter based on score. */
				if(0 == FilterAlignEntry(&a->entriesTwo[j],
							INT_MIN,
							startContig,
							startPos,
							endContig,
							endPos) &&
						(minScoreReadsPaired <= (a->entriesOne[i].score + a->entriesTwo[j].score))) {
					/* Get the distance between the two */
					curContigDistance = a->entriesTwo[j].contig - a->entriesOne[i].contig;
					curPositionDistance = a->entriesTwo[j].position - a->entriesOne[i].position;
					/* Get the distance from the mean */
					curDistanceToMean = abs(curPositionDistance - meanDistancePaired);
					if(0==curContigDistance &&
							curPositionDistance <= maxDistancePaired &&
							curPositionDistance >= minDistancePaired) {
						/* Inside the bounds */

						switch(algorithmReadsPaired) {
							case Unique:
								/* unique */
								numUnique++;
								uniqueOne = i;
								uniqueTwo = j;
								break;
							case BestScore:
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
								break;
							case MeanUnique:
								/* closest unique */
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
								break;
							case MeanBestScore:
								/* closest best score */
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
								break;
							default:
								PrintError(FnName,
										"algorithmReadsPaired",
										"Could not understand algorithmReadsPaired",
										Exit,
										OutOfRange);
								break;
						}
					}
					else {
						/* Outside the bounds */

						switch(algorithmReadsPaired) {
							case Unique:
							case MeanUnique:
								/* unique outside */
								numUniqueOutside++;
								uniqueOutsideOne = i;
								uniqueOutsideTwo = j;
								break;
							case BestScore:
							case MeanBestScore:
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
								break;
							default:
								PrintError(FnName,
										"algorithmReadsPaired",
										"Could not understand algorithmReadsPaired",
										Exit,
										OutOfRange);
								break;
						}
					}
				}
			}
		}
	}

	/* Check if we have a unique pair of alignments */
	switch(algorithmReadsPaired) {
		case Unique: /* unique */
			if(numUnique == 1) { 
				/* Prefer within the boundaries */
				AlignEntriesKeepOnly(a,
						uniqueOne,
						uniqueTwo,
						1,
						a->space);
				foundType=Found;
			}
			else if(numUnique == 0 &&
					numUniqueOutside == 1) { 
				/* If there are no pairs within the bounds */
				AlignEntriesKeepOnly(a,
						uniqueOutsideOne,
						uniqueOutsideTwo,
						1,
						a->space);
				foundType=OutsideBounds;
			}
			else {
				foundType=NoneFound;
			}
			break;
		case BestScore: /* best score */
			if(numBestScore == 1) {
				/* Prefer within the boundaries */
				AlignEntriesKeepOnly(a,
						bestScoreOne,
						bestScoreTwo,
						1,
						a->space);
				foundType=Found;
			}
			else if(numBestScore == 0 &&
					numBestScoreOutside == 1) {
				/* If there are no pairs within the bounds */
				AlignEntriesKeepOnly(a,
						bestScoreOutsideOne,
						bestScoreOutsideTwo,
						1,
						a->space);
				foundType=OutsideBounds;
			}
			else {
				foundType=NoneFound;
			}
			break;
		case MeanUnique: /* closest to the mean unique */
		case MeanBestScore: /* closest to the mean best score */
			if(numClosestToMean == 1) {
				AlignEntriesKeepOnly(a,
						closestToMeanOne,
						closestToMeanTwo,
						1,
						a->space);
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
			foundType = ContigAb;
		}
		else {
			/* Different strand - inversion */
			foundType = Inversion;
		}
	}
	return foundType;
}
