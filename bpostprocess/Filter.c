#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <limits.h>
#include <float.h>
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
		int maxMismatches,
		int maxColorErrors,
		int pairedEnd,
		int algorithmReadsPaired,
		int minScoreReadsPaired,
		int minDistancePaired,
		int maxDistancePaired,
		int maxMismatchesPaired,
		int maxColorErrorsPaired)
{
	char *FnName="FilterAlignEntries";
	int foundType=NoneFound;
	AlignEntries tmpA;
	assert(a->pairedEnd == pairedEnd);
	AlignEntriesInitialize(&tmpA);

	/* We should only modify "a" if it is going to be reported */ 

	/* Copy in case we do not find anything to report */
	AlignEntriesCopy(a, &tmpA);

	/* Filter for reads that are not paired end */
	if(pairedEnd == SingleEnd) {
		switch(algorithmReads) {
			case NoFiltering:
				/* Do nothing */
				foundType=Found;
				break;
			case AllNotFiltered:
				foundType=FilterReadInAlignEntries(a,
						algorithmReads,
						minScoreReads,
						startContig,
						startPos,
						endContig,
						endPos,
						maxMismatches,
						maxColorErrors,
						First);
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
						maxMismatches,
						maxColorErrors,
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
				if(Found==FilterReadInAlignEntries(a,
							algorithmReadsPaired,
							minScoreReads,
							startContig,
							startPos,
							endContig,
							endPos,
							maxMismatches,
							maxColorErrors,
							First) &&
						Found==FilterReadInAlignEntries(a,
							algorithmReadsPaired,
							minScoreReads,
							startContig,
							startPos,
							endContig,
							endPos,
							maxMismatches,
							maxColorErrors,
							Second)) {
					foundType=Found;
				}
				break;
			case Unique:
			case BestScore:
				foundType=FilterOneAlignEntries(a,
						startContig,
						startPos,
						endContig,
						endPos,
						minScoreReads,
						maxMismatches,
						maxColorErrors,
						algorithmReadsPaired,
						minScoreReadsPaired,
						minDistancePaired,
						maxDistancePaired,
						maxMismatchesPaired,
						maxColorErrorsPaired);
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
	/* If we filtered all, copy back */
	if(foundType == NoneFound) {
		AlignEntriesFree(a);
		AlignEntriesCopy(&tmpA, a);
	}
	AlignEntriesFree(&tmpA);

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
		int maxMismatches,
		int maxColorErrors,
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
					a->space,
					minScoreReads,
					startContig,
					startPos,
					endContig,
					endPos,
					maxMismatches,
					maxColorErrors,
					NULL,
					NULL)) {
			uniqueIndex=i;
			numNotFiltered++;

			switch(algorithmReads) {
				case AllNotFiltered:
					/* Do nothing */
					break;
				case Unique:
					if(numNotFiltered > 1) {
						/* This means that two ore more entries passed the filters.  Therefore
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
		int space,
		int minScoreReads,
		int startContig,
		int startPos,
		int endContig,
		int endPos,
		int maxMismatches,
		int maxColorErrors,
		int *returnNumMismatches,
		int *returnNumColorErrors)
{
	int32_t i, numMismatches, numColorErrors;
	numMismatches=numColorErrors=0;
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
	/* Check mismatches */
	for(i=0;i<a->length;i++) {
		/* Increment mismatches */
		if(GAP != a->read[i] && GAP != a->reference[i] && a->read[i] != a->reference[i]) {
			numMismatches++;
		}
	}
	if(NULL != returnNumMismatches) {
		(*returnNumMismatches) = numMismatches;
	}
	/* Check color errors */
	if(ColorSpace == space) {
		for(i=0;i<a->length;i++) {
			/* Increment color errors */
			switch(a->colorError[i]) {
				case '1':
					numColorErrors++;
					break;
				default:
					break;
			}
		}
	}
	if(NULL != returnNumColorErrors) {
		(*returnNumColorErrors) = numColorErrors;
	}
	if(maxMismatches < numMismatches) {
		return 1;
	}
	if(ColorSpace == space) {
		if(maxColorErrors < numColorErrors) {
			return 1;
		}
	}
	return 0;
}

/* TODO */
int FilterOneAlignEntries(AlignEntries *a,
		int startContig,
		int startPos,
		int endContig,
		int endPos,
		int minScoreReads,
		int maxMismatches,
		int maxColorErrors,
		int algorithmReadsPaired,
		int minScoreReadsPaired,
		int minDistancePaired,
		int maxDistancePaired,
		int maxMismatchesPaired,
		int maxColorErrorsPaired)
{
	char *FnName = "FilterOneAlignEntries";
	int numMismatchesOne, numMismatchesTwo, numColorErrorsOne, numColorErrorsTwo;
	/* best score */
	double bestScore=DBL_MIN;
	int bestScoreOne=-1, bestScoreTwo=-1;
	int numBestScore=0;
	/* unique */
	int uniqueOne=-1, uniqueTwo=-1;
	int numUnique=0;
	/* best score outside */
	double bestScoreOutside=DBL_MIN;
	int bestScoreOutsideOne=-1, bestScoreOutsideTwo=-1;
	int numBestScoreOutside=0;
	/* unique outside */
	int uniqueOutsideOne=-1, uniqueOutsideTwo=-1;
	int numUniqueOutside=0;
	/* current distance to the mean */
	int curContigDistance, curPositionDistance;
	double curScore;
	int i, j;
	int foundType=NoneFound;

	/* Note: this does not take into account strandedness */

	/* Go through all possible pairs of alignments */
	for(i=0;i<a->numEntriesOne;i++) { /* First read */
		/* Filter first entry for contig position */
				numMismatchesOne=numColorErrorsOne=0;
		if(0 == FilterAlignEntry(&a->entriesOne[i],
					a->space,
					minScoreReads,
					startContig,
					startPos,
					endContig,
					endPos,
					maxMismatches,
					maxColorErrors,
					&numMismatchesOne,
					&numColorErrorsOne)) {
			for(j=0;j<a->numEntriesTwo;j++) { /* Second read */
				numMismatchesTwo=numColorErrorsTwo=0;
				/* Get the current score */
				curScore = a->entriesOne[i].score + a->entriesTwo[j].score;
				/* Filter second entry for contig position.  Also filter based on score. */
				if(0 == FilterAlignEntry(&a->entriesTwo[j],
							a->space,
							minScoreReads,
							startContig,
							startPos,
							endContig,
							endPos,
							maxMismatches,
							maxColorErrors,
							&numMismatchesOne,
							&numColorErrorsOne) &&
						(minScoreReadsPaired <= curScore) &&
						(numMismatchesOne + numMismatchesTwo <= maxMismatchesPaired) &&
						(numColorErrorsOne + numColorErrorsTwo <= maxColorErrorsPaired)) {
					/* Get the distance between the two */
					curContigDistance = a->entriesTwo[j].contig - a->entriesOne[i].contig;
					curPositionDistance = a->entriesTwo[j].position - a->entriesOne[i].position;
					if(0==curContigDistance &&
							curPositionDistance <= maxDistancePaired &&
							curPositionDistance >= minDistancePaired) {
						/* Inside the bounds */
						switch(algorithmReadsPaired) {
							case Unique:
								/* unique */
								numUnique++; /* Update the number of unique */
								uniqueOne = i;
								uniqueTwo = j;
								break;
							case BestScore:
								/* best score */
								if(bestScore == curScore) {
									/* If we have the same score update the number */ 
									assert(numBestScore > 0);
									numBestScore++;
								}
								else if(bestScore < curScore) {
									/* Update the score */
									bestScore = curScore;
									bestScoreOne = i;
									bestScoreTwo = j;
									numBestScore=1;
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
								/* unique outside */
								numUniqueOutside++;
								uniqueOutsideOne = i;
								uniqueOutsideTwo = j;
								break;
							case BestScore:
								/* best score outside */
								if(curScore == bestScoreOutside) {
									/* If we have the same score update the number */ 
									assert(numBestScoreOutside > 0);
									numBestScoreOutside++;
								}
								else if(bestScoreOutside < curScore) {
									/* Update the score */
									bestScoreOutside = curScore;
									bestScoreOutsideOne = i;
									bestScoreOutsideTwo = j;
									numBestScoreOutside=1;
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
				/* Free */
				AlignEntriesFree(a);
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
				/* Free */
				AlignEntriesFree(a);
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
			/* Same strand - contig abnormality */
			foundType = ContigAb;
		}
		else {
			/* Different strand - inversion */
			foundType = Inversion;
		}
	}
	return foundType;
}
