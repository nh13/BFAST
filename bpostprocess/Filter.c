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
		int unpaired,
		int minScoreReadsPaired,
		int minDistancePaired,
		int maxDistancePaired,
		int maxMismatchesPaired,
		int maxColorErrorsPaired)
{
	char *FnName="FilterAlignEntries";
	int foundType=NoneFound;
	AlignEntries tmpA;
	int32_t i, bestScore, bestScoreIndex, numBestScore;
	assert(a->pairedEnd == pairedEnd);

	AlignEntriesInitialize(&tmpA);

	if(SingleEnd == a->pairedEnd && 
			NoFiltering == algorithmReads) {
		return Found;
	}
	else if(PairedEnd == a->pairedEnd &&
			NoFiltering == algorithmReadsPaired) {
		return Found;
	}

	/* We should only modify "a" if it is going to be reported */ 
	/* Copy in case we do not find anything to report */
	AlignEntriesCopy(a, &tmpA);

	/* Filter each read individually */
	for(i=0;i<tmpA.numEntriesOne;i++) {
		/* Check if we should filter */
		if(0<FilterAlignEntry(&tmpA.entriesOne[i],
					tmpA.space,
					minScoreReads,
					startContig,
					startPos,
					endContig,
					endPos,
					maxMismatches,
					maxColorErrors)) {
			/* Copy end here */
			if(i < tmpA.numEntriesOne-1) {
				AlignEntryCopy(&tmpA.entriesOne[tmpA.numEntriesOne-1], &tmpA.entriesOne[i]);
			}
			AlignEntriesReallocate(&tmpA, tmpA.numEntriesOne-1, tmpA.numEntriesTwo, tmpA.pairedEnd, tmpA.space);
			/* Since we removed, do not incrememt */
			i--;
		}
	}
	if(PairedEnd == tmpA.pairedEnd) {
		for(i=0;i<tmpA.numEntriesTwo;i++) {
			/* Check if we should filter */
			if(0<FilterAlignEntry(&tmpA.entriesTwo[i],
						tmpA.space,
						minScoreReads,
						startContig,
						startPos,
						endContig,
						endPos,
						maxMismatches,
						maxColorErrors)) {
				/* Copy end here */
				if(i < tmpA.numEntriesTwo-1) {
					AlignEntryCopy(&tmpA.entriesTwo[tmpA.numEntriesTwo-1], &tmpA.entriesTwo[i]);
				}
				AlignEntriesReallocate(&tmpA, tmpA.numEntriesOne, tmpA.numEntriesTwo-1, tmpA.pairedEnd, tmpA.space);
				/* Since we removed, do not incrememt */
				i--;
			}
		}
	}

	/* Pick alignment if possible */
	if(pairedEnd == SingleEnd) {
		switch(algorithmReads) {
			case AllNotFiltered:
				foundType=(0<tmpA.numEntriesOne)?Found:NoneFound;
				break;
			case Unique:
				foundType=(1==tmpA.numEntriesOne)?Found:NoneFound;
				break;
			case BestScore:
				bestScore = INT_MIN;
				bestScoreIndex = -1;
				numBestScore=0;
				for(i=0;i<tmpA.numEntriesOne;i++) {
					if(bestScore < tmpA.entriesOne[i].score) {
						bestScore = tmpA.entriesOne[i].score;
						bestScoreIndex = i;
						numBestScore = 1;
					}
					else if(bestScore == tmpA.entriesOne[i].score) {
						numBestScore++;
					}
				}
				if(1 == numBestScore) {
					AlignEntryCopy(&tmpA.entriesOne[bestScoreIndex], &tmpA.entriesOne[0]);
					AlignEntriesReallocate(a, 1, tmpA.numEntriesTwo, tmpA.pairedEnd, tmpA.space);
					foundType=Found;
				}
				else {
					foundType=NoneFound;
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
		if(tmpA.numEntriesOne == 0 && 0 < tmpA.numEntriesTwo && 1 == unpaired) {
			switch(algorithmReadsPaired) {
				case AllNotFiltered:
					foundType=(0<tmpA.numEntriesTwo)?Unpaired:NoneFound;
					break;
				case Unique:
					foundType=(1==tmpA.numEntriesTwo)?Unpaired:NoneFound;
					break;
				case BestScore:
					bestScore = INT_MIN;
					bestScoreIndex = -1;
					numBestScore=0;
					for(i=0;i<tmpA.numEntriesTwo;i++) {
						if(bestScore < tmpA.entriesTwo[i].score) {
							bestScore = tmpA.entriesTwo[i].score;
							bestScoreIndex = i;
							numBestScore = 1;
						}
						else if(bestScore == tmpA.entriesTwo[i].score) {
							numBestScore++;
						}
					}
					if(1 == numBestScore) {
						AlignEntryCopy(&tmpA.entriesTwo[bestScoreIndex], &tmpA.entriesTwo[0]);
						AlignEntriesReallocate(a, tmpA.numEntriesOne, 1, tmpA.pairedEnd, tmpA.space);
						foundType=Unpaired;
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
		}
		else if(0 < tmpA.numEntriesOne && tmpA.numEntriesTwo == 0 && 1 == unpaired) {
			switch(algorithmReadsPaired) {
				case AllNotFiltered:
					foundType=(0<tmpA.numEntriesOne)?Unpaired:NoneFound;
					break;
				case Unique:
					foundType=(1==tmpA.numEntriesOne)?Unpaired:NoneFound;
					break;
				case BestScore:
					bestScore = INT_MIN;
					bestScoreIndex = -1;
					numBestScore=0;
					for(i=0;i<tmpA.numEntriesOne;i++) {
						if(bestScore < tmpA.entriesOne[i].score) {
							bestScore = tmpA.entriesOne[i].score;
							bestScoreIndex = i;
							numBestScore = 1;
						}
						else if(bestScore == tmpA.entriesOne[i].score) {
							numBestScore++;
						}
					}
					if(1 == numBestScore) {
						AlignEntryCopy(&tmpA.entriesOne[bestScoreIndex], &tmpA.entriesOne[0]);
						AlignEntriesReallocate(a, 1, tmpA.numEntriesTwo, tmpA.pairedEnd, tmpA.space);
						foundType=Unpaired;
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
		}
		else {
			switch(algorithmReadsPaired) {
				case AllNotFiltered:
				case Unique:
				case BestScore:
					foundType = FilterPairedEnd(&tmpA,
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
	}
	/* If we found, then copy back */
	if(NoneFound != foundType) {
		AlignEntriesFree(a);
		AlignEntriesCopy(&tmpA, a);
	}
	AlignEntriesFree(&tmpA);

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
		int maxColorErrors)
{
	int32_t numMismatches, numColorErrors;
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
		return 2;
	}
	/* Check mismatches */
	numMismatches = GetNumMismatchesInAlignEntry(a);
	/* Check color errors */
	numColorErrors = GetNumColorErrorsInAlignEntry(a, space);
	if(maxMismatches < numMismatches) {
		return 3;
	}
	if(ColorSpace == space) {
		if(maxColorErrors < numColorErrors) {
			return 4;
		}
	}
	return 0;
}

/* TODO */
int GetNumMismatchesInAlignEntry(AlignEntry *a)
{
	int32_t i, numMismatches;
	/* Check mismatches */
	numMismatches=0;
	for(i=0;i<a->length;i++) {
		/* Increment mismatches */
		if(GAP != a->read[i] && 
				GAP != a->reference[i] && 
				a->read[i] != a->reference[i]) {
			numMismatches++;
		}
	}
	return numMismatches;
}

/* TODO */
int GetNumColorErrorsInAlignEntry(AlignEntry *a, int space)
{
	int32_t i, numColorErrors;
	/* Check color errors */
	numColorErrors=0;
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
	return numColorErrors;
}

/* TODO */
int FilterPairedEnd(AlignEntries *a,
		int algorithmReadsPaired,
		int minScoreReadsPaired,
		int minDistancePaired,
		int maxDistancePaired,
		int maxMismatchesPaired,
		int maxColorErrorsPaired)
{
	char *FnName = "FilterPairedEnd";
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
	int i, j, foundType;

	foundType=NoneFound;

	/* Go through all possible pairs of alignments */
	for(i=0;i<a->numEntriesOne;i++) { /* First read */
		/* Filter first entry for contig position */
		numMismatchesOne=GetNumMismatchesInAlignEntry(&a->entriesOne[i]);
		numColorErrorsOne=GetNumColorErrorsInAlignEntry(&a->entriesOne[i], a->space);
		for(j=0;j<a->numEntriesTwo;j++) { /* Second read */
			numMismatchesTwo=GetNumMismatchesInAlignEntry(&a->entriesTwo[j]);
			numColorErrorsTwo=GetNumColorErrorsInAlignEntry(&a->entriesTwo[j], a->space);
			/* Get the current score */
			curScore = a->entriesOne[i].score + a->entriesTwo[j].score;
			/* Filter second entry for contig position.  Also filter based on score. */
			if((minScoreReadsPaired <= curScore) &&
					(numMismatchesOne + numMismatchesTwo <= maxMismatchesPaired) &&
					(numColorErrorsOne + numColorErrorsTwo <= maxColorErrorsPaired)) {
				/* Get the distance between the two */
				curContigDistance = a->entriesTwo[j].contig - a->entriesOne[i].contig;
				curPositionDistance = a->entriesTwo[j].position - a->entriesOne[i].position;
				/* Check within bounds */
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
