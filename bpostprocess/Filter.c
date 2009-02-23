#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <limits.h>
#include <float.h>
#include "../blib/BError.h"
#include "../blib/AlignEntry.h"
#include "../blib/AlignEntries.h"
#include "../blib/BLib.h"
#include "Definitions.h"
#include "Filter.h"

/* TODO */
int FilterAlignEntries(AlignEntries *a,
		int algorithmReads,
		int uniquenessScore,
		int minUniquenessScore,
		int minScoreReads,
		int startContig,
		int startPos,
		int endContig,
		int endPos,
		int maxMismatches,
		int maxColorErrors,
		int pairedEnd,
		int minDistancePaired,
		int maxDistancePaired,
		int useDistancePaired,
		int contigAbPaired,
		int inversionsPaired,
		int unpaired)
{
	char *FnName="FilterAlignEntries";
	int foundType, foundTypeOne, foundTypeTwo;
	AlignEntries tmpA, tmpB;
	int32_t i;
	int32_t bestScoreOne, bestScoreIndexOne, numBestScoreOne;
	int32_t bestScoreTwo, bestScoreIndexTwo, numBestScoreTwo;
	int32_t curContigDistance, curPositionDistance;
	double uniquenessScoreOne, uniquenessScoreTwo;
	assert(a->pairedEnd == pairedEnd);

	AlignEntriesInitialize(&tmpA);

	/* We should only modify "a" if it is going to be reported */ 
	/* Copy in case we do not find anything to report */
	AlignEntriesCopy(a, &tmpA);

	if(NoFiltering != algorithmReads) {
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
	}
	foundType=foundTypeOne=foundTypeTwo=NoneFound;

	/* Pick alignment */
	switch(algorithmReads) {
		case NoFiltering:
		case AllNotFiltered:
			foundTypeOne=(0<tmpA.numEntriesOne)?Found:NoneFound;
			foundTypeTwo=(0<tmpA.numEntriesTwo)?Found:NoneFound;
			if(1==uniquenessScore) {
				AlignEntriesCopy(a, &tmpB);
				for(i=0;i<tmpB.numEntriesOne;i++) {
					tmpA.entriesOne[i].score=GetUniquenessScore(tmpB.entriesOne,
							tmpB.numEntriesOne,
							i);
				}
				for(i=0;i<tmpB.numEntriesTwo;i++) {
					tmpA.entriesTwo[i].score=GetUniquenessScore(tmpB.entriesTwo,
							tmpB.numEntriesTwo,
							i);
				}
				AlignEntriesFree(&tmpB);
			}
			break;
		case Unique:
			foundTypeOne=(1==tmpA.numEntriesOne)?Found:NoneFound;
			foundTypeTwo=(1==tmpA.numEntriesTwo)?Found:NoneFound;
			if(1==uniquenessScore) {
				if(Found == foundTypeOne) {
					tmpA.entriesOne[0].score = 1;
				}
				if(Found == foundTypeTwo) {
					tmpA.entriesTwo[0].score = 1;
				}
			}
			break;
		case BestScore:
			bestScoreOne = bestScoreTwo = INT_MIN;
			bestScoreIndexOne = bestScoreIndexTwo = -1;
			numBestScoreOne = numBestScoreTwo = 0;
			/* One */
			for(i=0;i<tmpA.numEntriesOne;i++) {
				if(bestScoreOne < tmpA.entriesOne[i].score) {
					bestScoreOne = tmpA.entriesOne[i].score;
					bestScoreIndexOne = i;
					numBestScoreOne = 1;
				}
				else if(bestScoreOne == tmpA.entriesOne[i].score) {
					numBestScoreOne++;
				}
			}
			if(1==numBestScoreOne) {
				foundTypeOne = Found;
				if(1==uniquenessScore) {
					/* Calculate uniqueness score */
					uniquenessScoreOne=GetUniquenessScore(tmpA.entriesOne,
							tmpA.numEntriesOne,
							bestScoreIndexOne);
					if(uniquenessScoreOne < minUniquenessScore) {
						foundTypeOne = NoneFound;
					}
					else {
						tmpA.entriesOne[bestScoreIndexOne].score = uniquenessScoreOne;
					}
				}
				if(Found == foundTypeOne) {
					AlignEntryCopy(&tmpA.entriesOne[bestScoreIndexOne], &tmpA.entriesOne[0]);
					AlignEntriesReallocate(&tmpA, 1, tmpA.numEntriesTwo, tmpA.pairedEnd, tmpA.space);
				}
			}
			/* Two */
			if(PairedEnd == tmpA.pairedEnd) {
				for(i=0;i<tmpA.numEntriesTwo;i++) {
					if(bestScoreTwo < tmpA.entriesTwo[i].score) {
						bestScoreTwo = tmpA.entriesTwo[i].score;
						bestScoreIndexTwo = i;
						numBestScoreTwo = 1;
					}
					else if(bestScoreTwo == tmpA.entriesTwo[i].score) {
						numBestScoreTwo++;
					}
				}
				if(1==numBestScoreTwo) {
					foundTypeTwo = Found;
					if(1==uniquenessScore) {
						/* Calculate uniqueness score */
						uniquenessScoreTwo=GetUniquenessScore(tmpA.entriesTwo,
								tmpA.numEntriesTwo,
								bestScoreIndexTwo);
						if(uniquenessScoreTwo < minUniquenessScore) {
							foundTypeTwo = NoneFound;
						}
						else {
							tmpA.entriesTwo[bestScoreIndexTwo].score = uniquenessScoreTwo;
						}
					}
					if(Found == foundTypeTwo) {
						AlignEntryCopy(&tmpA.entriesTwo[bestScoreIndexTwo], &tmpA.entriesTwo[0]);
						AlignEntriesReallocate(&tmpA, tmpA.numEntriesOne, 1, tmpA.pairedEnd, tmpA.space);
					}
				}
			}
			break;
		default:
			PrintError(FnName,
					"algorithmReadsPaired",
					"Could not understand algorithmReads",
					Exit,
					OutOfRange);
			break;
	}

	if(SingleEnd == tmpA.pairedEnd) {
		foundType=foundTypeOne;
	}
	else {
		foundType=(Found==foundTypeOne && Found==foundTypeTwo)?Found:NoneFound;

		if(Found == foundType) {
			if(1 == useDistancePaired) {
				curContigDistance = a->entriesTwo[0].contig - a->entriesOne[0].contig;
				curPositionDistance = a->entriesTwo[0].position - a->entriesOne[0].position;

				if(0 != curContigDistance ||
						curPositionDistance < minDistancePaired ||
						maxDistancePaired < curPositionDistance) {
					if(a->entriesOne[0].strand == a->entriesTwo[0].strand) {
						if(1 == contigAbPaired) {
							/* Same strand - contig abnormality */
							foundType = ContigAb;
						}
						else {
							foundType = NoneFound;
						}
					}
					else {
						if(1 == inversionsPaired) { 
							/* Different strand - inversion */
							foundType = Inversion;
						}
						else {
							foundType = NoneFound;
						}
					}
				}
			}
		}
		else if(1 == unpaired) {
			if(Found == foundTypeOne) {
				assert(NoneFound==foundTypeTwo);
				foundType = Unpaired;
				AlignEntriesReallocate(&tmpA, 1, 0, tmpA.pairedEnd, tmpA.space);
			}
			else if(Found == foundTypeTwo) {
				assert(NoneFound==foundTypeOne);
				foundType = Unpaired;
				AlignEntriesReallocate(&tmpA, 0, 1, tmpA.pairedEnd, tmpA.space);
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
				ToLower(a->read[i]) != ToLower(a->reference[i])) {
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
				case GAP:
					break;
				default:
					numColorErrors++;
					break;
			}
		}
	}
	return numColorErrors;
}

/* TODO */
double GetUniquenessScore(AlignEntry *a,
		int32_t numEntries,
		int32_t bestScoreIndex) 
{
	int32_t nextBestScore, numNextBestScore;
	int32_t i;

	nextBestScore = INT_MIN;
	numNextBestScore = 0;

	for(i=0;i<numEntries;i++) {
		if(i != bestScoreIndex) {
			if(nextBestScore < a[i].score) {
				nextBestScore = a[i].score;
				numNextBestScore = 1;
			}
			else if(nextBestScore == a[i].score) {
				numNextBestScore++;
			}
		}
	}

	if(0 < numNextBestScore && 0 < nextBestScore) {
		return (a[bestScoreIndex].score
				/ (a[bestScoreIndex].score + numNextBestScore*nextBestScore));
	}
	else {
		return 1.0;
	}
}
